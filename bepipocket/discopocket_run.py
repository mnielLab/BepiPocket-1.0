### IMPORTS ###
from pathlib import Path
import sys
import shutil
import subprocess
import numpy as np
import pickle
import pandas as pd
from string import ascii_uppercase

from chai_lab.chai1 import run_inference

MODULE_DIR = Path(__file__).resolve().parent
sys.path.append(str(MODULE_DIR))
from fasta_utilities import read_accs_and_sequences_from_fasta
from general_functions import _run_complete, _wipe_dir
from biopdb_utilities import is_pdb_file, is_cif_file, cif_to_pdb

from restraint_utilities import abag_make_pocket_restraints, abag_lightpocket_hcdr3_restraints
from anarci_utilities import get_hcdr3_center_residue

DISCOTOPE3_MODELS = MODULE_DIR / "models" / "discotope3_models"


def run_discotope3_pdb(structure_file, discotope_outdir, outdir, multichain_mode=True, cpu_only=False):
    """
    Run DiscoTope-3.0 on a single structure file (PDB or CIF).

    Parameters
    ----------
    structure_file : str or Path
        Input structure path.
    discotope_outdir : str or Path
        Temporary DiscoTope output directory.
    outdir : str or Path
        Final directory where CSV outputs are copied.
    multichain_mode : bool
        Whether to use DiscoTope multichain mode.
    cpu_only : bool
        Whether to force DiscoTope CPU mode.

    Returns
    -------
    list[Path]
        Copied DiscoTope CSV output files in `outdir`.
    """
    structure_file = Path(structure_file)
    discotope_outdir = Path(discotope_outdir)
    outdir = Path(outdir)

    if not DISCOTOPE3_MODELS.is_dir():
        raise FileNotFoundError(f"DiscoTope model directory not found: {DISCOTOPE3_MODELS}")

    cmd = [
        sys.executable,
        "-m",
        "discotope3.main",
        "--models_dir",
        str(DISCOTOPE3_MODELS),
    ]

    if is_pdb_file(structure_file):
        pdb_input = structure_file
    elif is_cif_file(structure_file):
        pdb_input = structure_file.parent / f"{structure_file.stem}.pdb"
        cif_to_pdb(structure_file, pdb_input)
    else:
        raise ValueError(f"Specified structure file was not a valid PDB or CIF file: {structure_file}")

    cmd.extend(["-f", str(pdb_input)])

    if multichain_mode:
        cmd.append("--multichain_mode")
    if cpu_only:
        cmd.append("--cpu_only")

    cmd.extend(["--out_dir", str(discotope_outdir)])

    discotope_outdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    subprocess.run(cmd, check=True)

    discotope_csv_files = list(discotope_outdir.glob("**/*.csv"))
    copied_csvs = []

    for f in discotope_csv_files:
        dest = outdir / f.name
        shutil.copy(f, dest)
        copied_csvs.append(dest)

    if multichain_mode:
        discotope_done_file = outdir / "discotope3_multimer_done.txt"
    else:
        discotope_done_file = outdir / "discotope3_monomer_done.txt"
    discotope_done_file.touch()
    
    # clean up
    pdb_input.unlink()
    if discotope_outdir.exists():
        shutil.rmtree(discotope_outdir)

    return copied_csvs


def parse_discotope3_csv(csv_file):
    csv_file = Path(csv_file)
    df = pd.read_csv(csv_file)

    required_cols = {"chain", "res_id", "residue", "DiscoTope-3.0_score"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required DiscoTope columns in {csv_file}: {sorted(missing_cols)}")

    df["chain_id"] = df["chain"]

    # Sort by chain + residue index
    df = df.sort_values(["chain_id", "res_id"])

    sequences = []
    scores = []
    chain_ids = []

    for chain_id, group in df.groupby("chain_id", sort=False):
        sequences.append("".join(group["residue"].tolist()))
        scores.append(group["DiscoTope-3.0_score"].to_numpy(dtype=float))
        chain_ids.append(chain_id)

    return sequences, scores, chain_ids

def discopocket_run(
    fasta_path,
    outdir,
    num_trunk_recycles=4,
    num_diffn_timesteps=200,
    overwrite_earlier_jobcontent=False,
    nr_runs=5,
    max_distance_angstrom="10.0",
    msa_directory=None,
    num_ab_chains=2,
    hcdr3_mode=False,
    cpu_only_discotope=False,
):
    """
    DiscoPocket:
      1. Run Chai-1 on full complex
      2. Run Chai-1 on antigen-only
      3. Run DiscoTope-3.0 on best antigen-only structure
      4. Use DiscoTope-ranked residues to define restraint-guided reruns
    """

    outdir = Path(outdir)
    fasta_path = Path(fasta_path)
    outdir.mkdir(parents=True, exist_ok=True)

    #
    # 1. Initial Chai full-complex run
    #
    out_path = outdir / "seed0"
    if not _run_complete(out_path):

        if out_path.is_dir():
            _wipe_dir(out_path)

        run_inference(
            fasta_file=fasta_path,
            output_dir=out_path,
            num_trunk_recycles=num_trunk_recycles,
            num_diffn_timesteps=num_diffn_timesteps,
            seed=0,
            use_esm_embeddings=True,
            msa_directory=msa_directory,
        )

        (outdir / "initialrun_done.txt").touch()

    #
    # 2. Metadata / chain bookkeeping
    #
    accs_and_seqs = read_accs_and_sequences_from_fasta(fasta_path)
    
    nr_chains = len(accs_and_seqs)
   
    structure_letters = [ascii_uppercase[i] for i in range(nr_chains)]
    #
    antigen_letters = structure_letters[:-num_ab_chains]
    antigen_seqs = [d[1] for d in accs_and_seqs[:-num_ab_chains]]
    #
    antibody_letters = structure_letters[-num_ab_chains:]
    
    if hcdr3_mode:
        anarci_antigen_letters, light_chain_letter, heavy_chain_letter, center_hcdr3_residue = get_hcdr3_center_residue(
            fasta_path,
            outdir / "anarci_run",
            structure_letters,
        )

    # 4. Antigen-only Chai run
    antigen_only_fasta = outdir / "antigen_only.fasta"
    fasta_content = "\n".join(f">{acc}\n{seq}" for acc, seq in accs_and_seqs[:-num_ab_chains])
    with open(antigen_only_fasta, "w") as outfile: outfile.write(fasta_content)

    antigen_only_outdir = outdir / "antigen_only"
    if not _run_complete(antigen_only_outdir):
        if antigen_only_outdir.is_dir():
            _wipe_dir(antigen_only_outdir)

        run_inference(
            fasta_file=antigen_only_fasta,
            output_dir=antigen_only_outdir,
            num_trunk_recycles=num_trunk_recycles,
            num_diffn_timesteps=num_diffn_timesteps,
            seed=0,
            use_esm_embeddings=True,
            msa_directory=msa_directory,
        )

        (antigen_only_outdir / "done.txt").touch()

    antigen_score_files = list(antigen_only_outdir.glob("scores*"))
    antigen_scores = [np.load(score_file) for score_file in antigen_score_files]
    antigen_scores = np.asarray([0.8 * score["iptm"][0] + 0.2 * score["ptm"][0] for score in antigen_scores])

    best_ag_idx = np.argmax(antigen_scores)
    best_ag_scorefile = antigen_score_files[best_ag_idx]
    best_ag_structure = best_ag_scorefile.parent / f"{best_ag_scorefile.stem.replace('scores', 'pred')}.cif"

    #
    # 5. Run DiscoTope on antigen-only structure
    #
    discotope_run_outdir = outdir / "discotope3_tmp"
    discotope_final_outdir = outdir / "discotope3"

    discotope_csvs = run_discotope3_pdb(
        best_ag_structure,
        discotope_run_outdir,
        discotope_final_outdir,
        multichain_mode=(len(antigen_letters) > 1),
        cpu_only=cpu_only_discotope,
    )

    if len(discotope_csvs) == 0:
        raise RuntimeError("DiscoTope-3.0 did not produce any CSV outputs.")

    # should only be one .csv, but just in case
    if len(discotope_csvs) != 1:
        raise RuntimeError(
            f"Expected exactly 1 DiscoTope CSV output, found {len(discotope_csvs)}: {discotope_csvs}"
        )

    discotope_csv = discotope_csvs[0]
    ag_disco3_seqs, disco3_scores, disco3_chain_ids = parse_discotope3_csv(discotope_csv)

    if antigen_seqs != ag_disco3_seqs:
        raise ValueError(
            f"Antigen sequences from FASTA and DiscoTope do not match:\n"
            f"FASTA: {antigen_seqs}\n"
            f"DiscoTope: {ag_disco3_seqs}"
        )

    if antigen_letters != disco3_chain_ids:
        raise ValueError(
            f"Antigen chain IDs from FASTA and DiscoTope do not match:\n"
            f"FASTA: {antigen_letters}\n"
            f"DiscoTope: {disco3_chain_ids}"
        )

    ag_residue_score_lookup = {}
    for seq, chain_id, chain_scores in zip(ag_disco3_seqs, disco3_chain_ids, disco3_scores):
        if len(seq) != len(chain_scores):
            raise ValueError(
                f"Length mismatch for chain {chain_id}: "
                f"{len(seq)} residues but {len(chain_scores)} DiscoTope scores"
            )

        for i, (aa, score) in enumerate(zip(seq, chain_scores)):
            ag_residue_score_lookup[(aa, chain_id, i)] = float(score)
    
    discotope_score_lookup_path = outdir / "discotope_score_lookup.pickle"
    with open(discotope_score_lookup_path, "wb") as outfile:
        pickle.dump(ag_residue_score_lookup, outfile)

    sorted_antigen_residue_list = [
        k for k, _ in sorted(ag_residue_score_lookup.items(), key=lambda item: item[1], reverse=True)
    ]

    # 7. Iterative restraint-guided Chai runs
    restraintsdir = outdir / "restraints"
    restraintsdir.mkdir(parents=True, exist_ok=True)

    max_runs = min(nr_runs - 1, len(sorted_antigen_residue_list))
    for i in range(max_runs):

        out_path = outdir / f"discopocket{i}"

        if out_path.is_dir() and overwrite_earlier_jobcontent:
            _wipe_dir(out_path)

        run_file_check = _run_complete(out_path)

        if run_file_check:
            print(f"Skipping. Found all structure and confidence files for run {i} at {out_path}.")
            continue

        if out_path.is_dir() and not run_file_check:
            _wipe_dir(out_path)

        pred_epitope_residue = sorted_antigen_residue_list[i]
        pred_epitope_residues = [pred_epitope_residue]

        # sequence index -> restraint index (1-based)
        pred_epitope_residues = [(e[0], e[1], e[2] + 1) for e in pred_epitope_residues]

        restraint_file = restraintsdir / f"discopocket{i}.restraints"

        if hcdr3_mode:
            abag_lightpocket_hcdr3_restraints(
                pred_epitope_residues,
                restraint_file,
                light_chain_letter,
                center_hcdr3_residue,
                confidence="1.0",
                min_distance_angstrom="0.0",
                max_distance_angstrom="10.0",
            )
        else:
            abag_make_pocket_restraints(
                pred_epitope_residues,
                restraint_file,
                antibody_letters,
                max_distance_angstrom=max_distance_angstrom,
            )

        run_inference(
            fasta_file=fasta_path,
            output_dir=out_path,
            constraint_path=restraint_file,
            num_trunk_recycles=num_trunk_recycles,
            num_diffn_timesteps=num_diffn_timesteps,
            seed=0,
            use_esm_embeddings=True,
            msa_directory=msa_directory,
        )

    (outdir / "done.txt").touch()
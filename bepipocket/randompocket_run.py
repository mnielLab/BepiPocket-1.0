### IMPORTS ###

from pathlib import Path
import sys
MODULE_DIR = str(Path(Path(__file__).parent.resolve()))
sys.path.append(MODULE_DIR)
import numpy as np
import string
import pickle
#from chai_lab.chai1 import run_inference
from fasta_utilities import read_accs_and_sequences_from_fasta
from biopdb_utilities import get_epitope_patch_residues, collect_epitope_contacts
from general_functions import load_pickle_file, _run_complete, _wipe_dir
from restraint_utilities import abag_make_pocket_restraints

### STATIC VARIABLES ###

ascii_uppercase = list(string.ascii_uppercase)

### FUNCTIONS ###

def run_randomscore_fasta(fasta_path, seed=0):
    """
    Generate reproducible random residue scores in [0, 1) for each sequence in the FASTA.
    Returns:
        sequences: list of sequences
        random_scores: list of numpy arrays, one per sequence
    """
    accs_and_seqs = read_accs_and_sequences_from_fasta(fasta_path)
    sequences = [seq for _, seq in accs_and_seqs]

    rng = np.random.default_rng(seed)
    random_scores = [rng.random(len(seq)) for seq in sequences]

    return sequences, random_scores


def randompocket_run(
    fasta_path,
    outdir,
    random_score_lookup=None,
    random_seed=0,
    num_trunk_recycles=4,
    num_diffn_timesteps=200,
    overwrite_earlier_jobcontent=False,
    antigen_seqidxs=None,
    nr_runs=5,
    patch_mode=False,
    max_distance_angstrom="10.0",
    patch_angradius=6.0,
    max_patch_size=4,
    epipara_aang_distance=5,
    msa_directory=None,
    num_ab_chains=2,
):
    """
    fasta_path: Filepath to the fasta file. Both antigen and antibody chains are expected.
                Expects that the first entries are antigen chains.
                And the last entries are antibody chains.
                If there is both a light and heavy chains, the two last entries should be for these two chains.
                Their order does not matter.
                If there is only a single antibody chain (nanobody), the last entry needs to be this chain.

    num_ab_chains: Default is 2 (light and heavy chain). If working with nanobody, set this to '1'.

    num_trunk_recycles=4, num_diffn_timesteps=200 are the default settings from Chai-1
    runs=None means to run until entire surface accessibility area is covered.
    epipara_ang_distance=5Å (using the same Å distance used to define contacts for DockQ)

    random_seed: seed for reproducible random residue scores
    """

    outdir = Path(outdir)
    
    out_path = outdir / "seed0"
    if not _run_complete(out_path):

        # if run terminated prematurely, re-run
        if out_path.is_dir(): _wipe_dir(out_path)


        ## initial run with no surface area restraint ##
        run_inference(
            fasta_file=fasta_path,
            output_dir=out_path,
            num_trunk_recycles=num_trunk_recycles,
            num_diffn_timesteps=num_diffn_timesteps,
            seed=0,
            use_esm_embeddings=True,
            msa_directory=msa_directory,
        )

        # write initial run done file
        outfile = open(outdir / "initialrun_done.txt", "w")
        outfile.close()

    # get epitope contacts of initial structure
    init_epitope_contacts, _ = collect_epitope_contacts(
        list(Path(outdir / "seed0").glob("pred*")),
        epipara_aang_distance=epipara_aang_distance,
    )
    all_epitope_contacts = init_epitope_contacts

    # antibody and antigen letters
    accs_and_seqs = read_accs_and_sequences_from_fasta(fasta_path)
    nr_chains = len(accs_and_seqs)
    structure_letters = [ascii_uppercase[i] for i in range(nr_chains)]
    antigen_letters, antibody_letters = structure_letters[:-num_ab_chains], structure_letters[-num_ab_chains:]

    aa_chainletter_residxs = []
    for i in range(nr_chains):
        acc, seq = accs_and_seqs[i]
        structure_letter = structure_letters[i]
        aa_chainletter_residx = [(aa, structure_letter, j) for j, aa in enumerate(seq)]
        aa_chainletter_residxs.append(aa_chainletter_residx)

    # get antigen fasta entries
    if antigen_seqidxs is None:
        ag_accs_and_seqs = [acc_seq for acc_seq in accs_and_seqs[:-num_ab_chains]]
        ag_aa_chainletter_residxs = [x for x in aa_chainletter_residxs[:-num_ab_chains]]
    else:
        ag_accs_and_seqs = [accs_and_seqs[idx] for idx in antigen_seqidxs]
        ag_aa_chainletter_residxs = [aa_chainletter_residxs[idx] for idx in antigen_seqidxs]

    ag_seqs = [d[1] for d in ag_accs_and_seqs]

    # identify highest confidence structure
    out_path = outdir / "seed0"
    score_files = list(out_path.glob("scores*"))
    scores = [np.load(score_file) for score_file in score_files]
    scores = np.asarray([0.8 * score["iptm"][0] + 0.2 * score["ptm"][0] for score in scores])
    best_scoreidx = np.argmax(scores)
    best_score, best_scorefile = scores[best_scoreidx], score_files[best_scoreidx]
    rank1_structure = best_scorefile.parent / f"{best_scorefile.stem.replace('scores', 'pred')}.cif"
    print(f"Structure with highest confidence {rank1_structure} with {best_score}")

    # random score lookup from user input
    random_score_lookup_path = outdir / "random_score_lookup.pickle"
    if random_score_lookup is not None:
        random_score_lookup = load_pickle_file(random_score_lookup)
    # random score lookup from previous run
    elif random_score_lookup_path.is_file():
        random_score_lookup = load_pickle_file(random_score_lookup_path)
    # compute random scores
    else:
        print(f"Computing random scores with seed={random_seed}")
        scored_seqs, random_scores = run_randomscore_fasta(fasta_path, seed=random_seed)
        random_score_lookup = {k: v for k, v in zip(scored_seqs, random_scores)}

    # save random score run scores
    if not random_score_lookup_path.is_file():
        with open(random_score_lookup_path, "wb") as outfile:
            pickle.dump(random_score_lookup, outfile)

    # random scores not found, recompute
    scored_sequences = set(random_score_lookup.keys())
    if len(set(ag_seqs) - scored_sequences) != 0:
        print(f"For {fasta_path}, random scores for some antigen sequences were not found in lookup. Recomputing scores")
        scored_seqs, random_scores = run_randomscore_fasta(fasta_path, seed=random_seed)
        random_score_lookup = {k: v for k, v in zip(scored_seqs, random_scores)}
        with open(random_score_lookup_path, "wb") as outfile:
            pickle.dump(random_score_lookup, outfile)

    ag_aa_chainletter_residxs_random_scores = {}
    N1 = len(ag_seqs)
    for i in range(N1):
        ag_seq = ag_seqs[i]
        residue_scores = random_score_lookup[ag_seq]
        N2 = len(ag_seq)
        for j in range(N2):
            ag_aa_chainletter_residxs_random_scores[ag_aa_chainletter_residxs[i][j]] = residue_scores[j]

    # sort antigen residues by random score
    sorted_antigen_residue_list = [
        k[0]
        for k in sorted(
            ag_aa_chainletter_residxs_random_scores.items(),
            key=lambda item: item[1],
            reverse=True,
        )
    ]

    restraintsdir = outdir / "restraints"
    if not restraintsdir.is_dir():
        restraintsdir.mkdir(parents=True)

    max_runs = min(nr_runs - 1, len(sorted_antigen_residue_list))
    for i in range(max_runs):

        out_path = outdir / f"randommap{i}"

        # delete results from earlier run
        if out_path.is_dir() and overwrite_earlier_jobcontent:
            _wipe_dir(out_path)

        # check how many files are in seed
        run_file_check = _run_complete(out_path)
      
         # already completed run, skip
        if run_file_check:
            print(f"Skipping. Found all structure and conf. files for run {i} at {out_path}.")
            continue

        # if run terminated prematurely, re-run
        if out_path.is_dir() and not run_file_check: _wipe_dir(out_path)
        
        pred_epitope_residue = sorted_antigen_residue_list[i]
        pred_epitope_residues = [pred_epitope_residue]

        # adjust residue indexing (PDB index starts at 1)
        pred_epitope_residues = [(e[0], e[1], e[2] + 1) for e in pred_epitope_residues]

        restraint_file = restraintsdir / f"randommap{i}.restraints"
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

    outfile = open(outdir / "done.txt", "w")
    outfile.close()
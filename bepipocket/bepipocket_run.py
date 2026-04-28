### IMPORTS ###
from pathlib import Path
import sys
MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)
import string
import pdb
import pickle
import torch
from chai_lab.chai1 import run_inference
from fasta_utilities import read_accs_and_sequences_from_fasta
from bp3 import bepipred3
from biopdb_utilities import prepare_epitope_patch_search, get_epitope_patch_residues
from general_functions import load_pickle_file, _run_complete, _wipe_dir, get_highest_confidence_structure
from restraint_utilities import abag_make_pocket_restraints, abag_lightpocket_hcdr3_restraints
from anarci_utilities import get_hcdr3_center_residue

### STATIC VARIABLES ###

ascii_uppercase = list(string.ascii_uppercase)

### FUNCTIONS ###

def get_antigen_data_from_fasta(fasta_path, num_ab_chains):
    accs_and_seqs = read_accs_and_sequences_from_fasta(fasta_path)
    nr_chains = len(accs_and_seqs)
    structure_letters = [ascii_uppercase[i] for i in range(nr_chains)] # will crash for (nr_chains > 26), very unlikely
    antibody_letters = structure_letters[-num_ab_chains:]

    aa_chainletter_residxs = []
    for i in range(nr_chains):
        seq = accs_and_seqs[i][1]
        structure_letter = structure_letters[i]
        aa_chainletter_residx = [(aa, structure_letter, j) for j, aa in enumerate(seq)]
        aa_chainletter_residxs.append(aa_chainletter_residx)

    
    if num_ab_chains < 1 or num_ab_chains >= len(accs_and_seqs):
        raise ValueError(f"num_ab_chains={num_ab_chains} is invalid for {len(accs_and_seqs)} chains")

    ag_accs_and_seqs = accs_and_seqs[:-num_ab_chains]
    ag_aa_chainletter_residxs = aa_chainletter_residxs[:-num_ab_chains]
    ag_seqs = [d[1] for d in ag_accs_and_seqs]

    return structure_letters, antibody_letters, ag_aa_chainletter_residxs, ag_seqs

def get_ag_residue_bp3_scores(fasta_path, outdir, ag_seqs, ag_aa_chainletter_residxs, bp3_score_lookup=None):
    bp3_score_lookup_path = outdir / "bp3_score_lookup.pickle"

    if bp3_score_lookup is not None:
        bp3_score_lookup = load_pickle_file(bp3_score_lookup)
    elif bp3_score_lookup_path.is_file():
        bp3_score_lookup = load_pickle_file(bp3_score_lookup_path)
    else:
        print("Computing BepiPred-3.0 scores")
        bp3scored_seqs, bp3_scores = run_bepipred3_fasta(fasta_path, outdir)
        bp3_score_lookup = {k: v for k, v in zip(bp3scored_seqs, bp3_scores)}

    if not bp3_score_lookup_path.is_file():
        with open(bp3_score_lookup_path, "wb") as outfile:
            pickle.dump(bp3_score_lookup, outfile)

    if len(set(ag_seqs) - set(bp3_score_lookup.keys())) != 0:
        print(f"For {fasta_path}, the BepiPred-3.0 scores for some antigen sequences were not found in lookup. Recomputing scores")
        bp3scored_seqs, bp3_scores = run_bepipred3_fasta(fasta_path, outdir, esm2_model_path=None)
        bp3_score_lookup = {k: v for k, v in zip(bp3scored_seqs, bp3_scores)}
        with open(bp3_score_lookup_path, "wb") as outfile:
            pickle.dump(bp3_score_lookup, outfile)

    ag_aa_chainletter_residxs_bp3_scores = {}
    for i, ag_seq in enumerate(ag_seqs):
        bp3_scores = bp3_score_lookup[ag_seq]
        for j, bp3_score in enumerate(bp3_scores):
            ag_aa_chainletter_residxs_bp3_scores[ag_aa_chainletter_residxs[i][j]] = bp3_score

    return ag_aa_chainletter_residxs_bp3_scores

def run_bepipred3_fasta(fasta_path, outdir, esm2_model_path=None, rm_esm2_encodings=True):

    # run BepiPred-3.0 on fasta
    esm2_enc_savepath = outdir / "esm2_encodings"
    antigens = bepipred3.Antigens(fasta_path,   esm2_enc_savepath, add_seq_len=True, run_esm_model_local=esm2_model_path)
    bp3_antigens_predict = bepipred3.BP3EnsemblePredict(antigens)
    bp3_antigens_predict.run_bp3_ensemble()
    antigen_sequences = antigens.seqs 

    # collect and average probabilities
    ensemble_probs = antigens.ensemble_probs
    avg_ensemble_probs = []
    for ensemble_prob in ensemble_probs:
        avg_ensemble_probs.append(torch.mean(torch.stack(ensemble_prob, axis=1), axis=1).cpu().detach().numpy())

    # remove temporary directory + files
    if rm_esm2_encodings:
        for f in esm2_enc_savepath.glob("*"): f.unlink()        
        esm2_enc_savepath.rmdir()


    return antigen_sequences, avg_ensemble_probs


def spread_epitope_ranking(sorted_antigen_residue_list, epitope_patch_lookup, ag_aa_chainletter_residxs_bp3_scores):
    new_ranked_list = []
    remaining_residues = sorted_antigen_residue_list[:]

    while remaining_residues:
        blocked_residues = set()
        next_round = []

        for ag_res in remaining_residues:
            if ag_res in blocked_residues:
                next_round.append(ag_res)
                continue

            new_ranked_list.append(ag_res)
            blocked_residues.update(epitope_patch_lookup[ag_res])

        remaining_residues = sorted(
            next_round,
            key=lambda ag_res: ag_aa_chainletter_residxs_bp3_scores[ag_res],
            reverse=True,
        )

    return new_ranked_list


def bepipocket_run(fasta_path, outdir, bp3_score_lookup=None, num_trunk_recycles=4, num_diffn_timesteps=200, overwrite_earlier_jobcontent=False,
                   nr_runs=5, max_distance_angstrom="10.0", epipara_aang_distance=5, msa_directory=None, num_ab_chains=2, hcdr3_mode=False, hobohm_patchradius=None):

    """
    fasta_path: Filepath to the fasta file. Both antigen and antibody chains are expected.
                Expects that the first entries are antigen chains.
                And the last entries are antibody chains.
                If there is both a light and heavy chains, the two last entries should be for these two chains.
                Their order does not matter.
                If there is only a single antibody chain (nanobody), the last entry needs to be this chain. 

    num_ab_chains: Default is 2 (light and heavy chain). If working with nanobody, set this to '1'.

    num_trunk_recycles=4, num_diffn_timesteps=200 are the default settings from Chai-1
    runs= None, means to run until entire surface accessibiltiy area is covered.
    epipara_ang_distance=5Å (using the same Å distance used to define contacts for DockQ)
    
    """
    outdir.mkdir(parents=True, exist_ok=True)
    
    out_path = outdir / "seed0"
    if not _run_complete(out_path):

        # if run terminated prematurely, re-run
        if out_path.is_dir(): _wipe_dir(out_path)

        ## initial run with no surface area restraint ##     
        run_inference(fasta_file=fasta_path, output_dir=out_path,
                          num_trunk_recycles=num_trunk_recycles,
                          num_diffn_timesteps=num_diffn_timesteps,
                          seed=0, use_esm_embeddings=True, msa_directory=msa_directory) 

        # write initial run done file 
        (outdir / "initialrun_done.txt").write_text("")

    # get antigen data + antibody letters from fasta
    structure_letters, antibody_letters, ag_aa_chainletter_residxs, ag_seqs = get_antigen_data_from_fasta(fasta_path, num_ab_chains)
    # use anarci to center hcdr3 
    if hcdr3_mode: _, light_chain_letter, _, center_hcdr3_residue = get_hcdr3_center_residue(fasta_path, outdir / "anarci_run", structure_letters)

    # get BepiPred-3.0 sorted score list
    ag_aa_chainletter_residxs_bp3_scores = get_ag_residue_bp3_scores(fasta_path, outdir, ag_seqs, ag_aa_chainletter_residxs, bp3_score_lookup=bp3_score_lookup)
    sorted_antigen_residue_list = [k[0] for k in sorted(ag_aa_chainletter_residxs_bp3_scores.items(), key=lambda item: item[1], reverse=True)]
    
    # diversify epitope ranking 
    if hobohm_patchradius is not None:
        rank1_structure_path = get_highest_confidence_structure(out_path)
        residues_by_chain, residue_index_lookup, search_atoms = prepare_epitope_patch_search(rank1_structure_path, num_ab_chains=num_ab_chains)
        epitope_patch_lookup = {ag_res: get_epitope_patch_residues(residues_by_chain, residue_index_lookup, search_atoms, ag_res, patch_angradius=hobohm_patchradius) for ag_res in sorted_antigen_residue_list}
        sorted_antigen_residue_list = spread_epitope_ranking(sorted_antigen_residue_list, epitope_patch_lookup, ag_aa_chainletter_residxs_bp3_scores)
    
    # run bepipocket for nr_runs (including initial run)
    restraintsdir = outdir / "restraints"
    if not restraintsdir.is_dir(): restraintsdir.mkdir(parents=True)
    max_runs = min(nr_runs - 1, len(sorted_antigen_residue_list))
    for i in range(max_runs):

        out_path = outdir / f"bepipocket{i}"

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
            
        # using one antigen residue as restraint
        pred_epitope_residue = sorted_antigen_residue_list[i] 
        pred_epitope_residues = [pred_epitope_residue]
        
        # adjust residue indexing (PDB index starts index 1. )
        pred_epitope_residues =  [(e[0], e[1], e[2]+1) for e in pred_epitope_residues]  
        
        # restraint outfile
        restraint_file = restraintsdir / f"bepipocket{i}.restraints" 
      
        if hcdr3_mode:
            abag_lightpocket_hcdr3_restraints(pred_epitope_residues, restraint_file, light_chain_letter, center_hcdr3_residue, confidence="1.0",
                                    min_distance_angstrom="0.0", max_distance_angstrom="10.0")

        # normal mode
        else:
            abag_make_pocket_restraints(pred_epitope_residues, restraint_file, antibody_letters, max_distance_angstrom=max_distance_angstrom)

        run_inference(fasta_file=fasta_path, output_dir=out_path,
                        constraint_path=restraint_file,
                        num_trunk_recycles=num_trunk_recycles,
                        num_diffn_timesteps=num_diffn_timesteps,
                        seed=0, use_esm_embeddings=True,
                        msa_directory=msa_directory)
 
    (outdir / "done.txt").write_text("")
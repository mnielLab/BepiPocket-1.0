### IMPORTS ###

from pathlib import Path
import sys
MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)
import numpy as np
import string
import pickle
import torch
from chai_lab.chai1 import run_inference
from fasta_utilities import read_accs_and_sequences_from_fasta
from bp3 import bepipred3
from biopdb_utilities import get_epitope_patch_residues, collect_epitope_contacts
from general_functions import load_pickle_file
from restraint_utilities import abag_make_pocket_restraints


### STATIC VARIABLES ###

ascii_uppercase = list(string.ascii_uppercase)

### FUNCTIONS ###

def run_bepipred3_fasta(fasta_path, outdir, esm2_model_path=None, rm_esm2_encodings=True):

    # run BepiPred-3.0 on fasta
    esm2_enc_savepath = outdir / "esm2_encodings"
    antigens = bepipred3.Antigens(fasta_path,   esm2_enc_savepath, add_seq_len=True, run_esm_model_local=esm2_model_path)
    bp3_antigens_predict = bepipred3.BP3EnsemblePredict(antigens)
    bp3_antigens_predict.run_bp3_ensemble()
    antigen_sequences = antigens.seqs 

    # collect and average probabilities
    ensemble_probs = antigens.ensemble_probs
    num_of_seqs = len(ensemble_probs)
    avg_ensemble_probs = []
    for ensemble_prob in ensemble_probs:
        avg_ensemble_probs.append(torch.mean(torch.stack(ensemble_prob, axis=1), axis=1).cpu().detach().numpy())

    # remove temporary directory + files
    if rm_esm2_encodings:
        for f in esm2_enc_savepath.glob("*"): f.unlink()        
        esm2_enc_savepath.rmdir()


    return antigen_sequences, avg_ensemble_probs
    

def bepipocket_run(fasta_path, outdir, bp3_score_lookup=None, num_trunk_recycles=4, num_diffn_timesteps=200,
                        overwrite_earlier_jobcontent=False, antigen_seqidxs=None, nr_runs=5, patch_mode=False, max_distance_angstrom="10.0",
                        patch_angradius=6.0, max_patch_size=4, epipara_aang_distance=5, msa_directory=None, num_ab_chains=2):

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

    if not Path(outdir / "initialrun_done.txt").is_file():

        ## initial run with no surface area restraint ##
        out_path = outdir / "seed0"
        run_inference(fasta_file=fasta_path, output_dir=out_path,
                          num_trunk_recycles=num_trunk_recycles,
                          num_diffn_timesteps=num_diffn_timesteps,
                          seed=0, use_esm_embeddings=True, msa_directory=msa_directory) 

        # write initial run done file 
        outfile = open(outdir / "initialrun_done.txt", "w")
        outfile.close()

    # get epitope contacts of initial structure
    init_epitope_contacts, _ = collect_epitope_contacts(list(Path(outdir / "seed0").glob("pred*")), epipara_aang_distance=epipara_aang_distance)
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
        ag_accs_and_seqs = [acc_seq for acc_seq in accs_and_seqs[:-num_ab_chains]] # assume antibody chains are always the last entries
        ag_aa_chainletter_residxs = [aa_chainletter_residx for aa_chainletter_residx in aa_chainletter_residxs[:-num_ab_chains]]
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

    # bepipred3 score from user input
    bp3_score_lookup_path = outdir / "bp3_score_lookup.pickle"
    if bp3_score_lookup is not None:
        bp3_score_lookup = load_pickle_file(bp3_score_lookup)
    # bepipred3 score from previous run
    elif bp3_score_lookup_path.is_file():
        bp3_score_lookup = load_pickle_file(bp3_score_lookup_path)
    # compute bepipred3 scores
    else:
        print("Computing BepiPred-3.0 scores")
        bp3scored_seqs, bp3_scores = run_bepipred3_fasta(fasta_path, outdir) 
        bp3_score_lookup = {k:v for k, v in zip(bp3scored_seqs, bp3_scores)}

    # save bepipred3 run scores
    if not bp3_score_lookup_path.is_file():
        with open(bp3_score_lookup_path, "wb") as outfile: pickle.dump(bp3_score_lookup, outfile)
    
    # bepipred-3.0 scores not found, recompute 
    bp3_scored_sequences = set(bp3_score_lookup.keys())
    if len( set(ag_seqs) - bp3_scored_sequences ) != 0:
        print(f"For {fasta_path}, the BepiPred-3.0 scores for some antigen sequences were not found in lookup. Recomputing scores")
        bp3scored_seqs, bp3_scores = run_bepipred3_fasta(fasta_path, outdir, esm2_model_path=None) 
        bp3_score_lookup = {k:v for k, v in zip(bp3scored_seqs, bp3_scores)}
        with open(bp3_score_lookup_path, "wb") as outfile: pickle.dump(bp3_score_lookup, outfile)

    ag_aa_chainletter_residxs_bp3_scores = {}
    N1 = len(ag_seqs)
    for i in  range(N1):
        ag_seq = ag_seqs[i]
        bp3_scores = bp3_score_lookup[ag_seq]
        N2 = len(ag_seq)
        for j in range(N2):
            ag_aa_chainletter_residxs_bp3_scores[ ag_aa_chainletter_residxs[i][j] ] = bp3_scores[j]

    # sort antigen residues by BepiPred-3.0 score
    sorted_antigen_residue_list = [k[0] for k in sorted(ag_aa_chainletter_residxs_bp3_scores.items(), key=lambda item: item[1], reverse=True)]
    
  
    restraintsdir = outdir / "restraints"
    if not restraintsdir.is_dir(): restraintsdir.mkdir(parents=True)
    
    for i in range(nr_runs - 1):

        out_path = outdir / f"bepipredmap{i}"

        # delete results from earlier run (if wanting to re-run evaluation)
        if out_path.is_dir() and overwrite_earlier_jobcontent:
            for f in out_path.glob("*"): f.unlink()

        # check how many files are in seed (should be 10: 5 .cif (structure) + 5 .npz (confidence))
        if out_path.is_dir():
            nr_score_files = len( list(out_path.glob("*.npz")) )
            nr_structure_files = len( list(out_path.glob("*.cif")) )
            if nr_structure_files == 5 and nr_score_files == 5: 
                run_file_check = True

        else: run_file_check = False

        if run_file_check:
            print(f"Skipping. Found all structure and conf. files for {i} at {str(out_path)}.")
            continue

        pred_epitope_residue = sorted_antigen_residue_list[i] 
         
        if patch_mode:
            # get epitope residues within XÅ of predicted epitope residue center
            epitope_residue_patch_residues = get_epitope_patch_residues(rank1_structure, pred_epitope_residue, patch_angradius=patch_angradius)
            
            # adjust index to match bepipred-3.0 score lookup
            epitope_residue_patch_residues = [(e[0], e[1], e[2]-1) for e in epitope_residue_patch_residues]
            
            # remove center residue
            idx = epitope_residue_patch_residues.index(pred_epitope_residue)
            epitope_residue_patch_residues = epitope_residue_patch_residues[:idx] + epitope_residue_patch_residues[idx+1:] 

            # sort surface patch residues by highest bepipred3 score
            res_bp3_scores = [(res, ag_aa_chainletter_residxs_bp3_scores[res]) for res in epitope_residue_patch_residues]
            res_bp3_scores = sorted(res_bp3_scores, key=lambda item: item[1], reverse=True)
            
            # if Y surface patch residues detected at X Å distance, less than max_patch_size, use Y surface patch residues
            patch_size = min([max_patch_size, len(res_bp3_scores)])

            # surface patch residues
            pred_epitope_residues = [pred_epitope_residue] + [res_bp3_scores[j][0] for j in range(patch_size - 1)]


        # just use one antigen residue as restraint (same as in preprint)
        else: pred_epitope_residues = [pred_epitope_residue]

        # adjust residue indexing (PDB index starts index 1. )
        pred_epitope_residues =  [(e[0], e[1], e[2]+1) for e in pred_epitope_residues]        
        restraint_file = restraintsdir / f"bepipredmap{i}.restraints" 
        abag_make_pocket_restraints(pred_epitope_residues, restraint_file, antibody_letters,
                                    max_distance_angstrom=max_distance_angstrom)

        run_inference(fasta_file=fasta_path, output_dir=out_path,
                        constraint_path=restraint_file,
                        num_trunk_recycles=num_trunk_recycles,
                        num_diffn_timesteps=num_diffn_timesteps,
                        seed=0, use_esm_embeddings=True,
                        msa_directory=msa_directory)
 
    outfile = open(outdir / "done.txt", "w")
    outfile.close()
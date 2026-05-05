### IMPORTS ###
from pathlib import Path
import sys
MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)

### FUNCTIONS ###
def abag_make_pocket_restraints(epitope_residues, restraint_outfile, antibody_chain_letters, confidence="1.0",
                                min_distance_angstrom="0.0", max_distance_angstrom="10.0", connection_type = "pocket"):
    """
    Generate pocket restraints for antibody-antigen interactions.

    Args:
        epitope_residues (list): List of epitope residues as tuples (res, chain_letter, residx).
        restraint_outfile (Path): Path to save the restraint file.
        antibody_chain_letters (list): List of antibody chain letters.
        confidence (str): Confidence level. Default is "1.0".
        min_distance_angstrom (str): Minimum distance in angstroms. Default is "0.0".
        max_distance_angstrom (str): Maximum distance in angstroms. Default is "10.0".
    """

    chai_restraint_header = "chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,max_distance_angstrom,comment,restraint_id"

    # Collect antigen restraints
    restraints = []
    for j, epi in enumerate(epitope_residues):

        res, antigen_chain_letter, residx = epi
        antigen_residx = f"{res}{residx}"

        # antibody chain restraints 
        for i, ab_chain in enumerate(antibody_chain_letters):
            ab_restraint = ",".join([
                ab_chain, "", antigen_chain_letter, antigen_residx, connection_type,
                confidence, min_distance_angstrom, max_distance_angstrom, f"antibodyc{i+1}-antigen"
            ])
            restraints.append(ab_restraint)

    # Add restraint ids
    nr_restraints = len(restraints)
    restraints = [f"{restraints[j]},restraint_{j}" for j in range(nr_restraints)]

    # Write restraint output file
    restraint_out = "\n".join([chai_restraint_header] + restraints)
    with open(restraint_outfile, "w") as outfile:
        outfile.write(restraint_out)


def abag_lightpocket_hcdr3_restraints(epitope_residues, restraint_outfile, light_chain_letter, hcdr3_residue, confidence="1.0",
                                min_distance_angstrom="0.0", max_distance_angstrom="10.0"):
    """
    Generate light-chain pocket restraints and HCDR3 contact restraints for antibody-antigen interactions.

    Args:
        epitope_residues (list): List of epitope residues as tuples (res, chain_letter, residx).
        restraint_outfile (Path): Path to save the restraint file.
        light_chain_letter (str or None): Light chain letter.
        hcdr3_residue (tuple): HCDR3 center residue as tuple (res, chain_letter, residx).
        confidence (str): Confidence level. Default is "1.0".
        min_distance_angstrom (str): Minimum distance in angstroms. Default is "0.0".
        max_distance_angstrom (str): Maximum distance in angstroms. Default is "10.0".
    """

    chai_restraint_header = "chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,max_distance_angstrom,comment,restraint_id"
    res, heavy_chain_letter, residx = hcdr3_residue 
    hcdr3_residx = f"{res}{residx}"

    # Collect antigen restraints
    restraints = []
    for j, epi in enumerate(epitope_residues):

        res, antigen_chain_letter, residx = epi
        antigen_residx = f"{res}{residx}"
        if light_chain_letter is not None:
            lightpocket_restraint = ",".join([
                    light_chain_letter, "", antigen_chain_letter, antigen_residx, "pocket",
                    confidence, min_distance_angstrom, max_distance_angstrom, f"antibodycL-antigen"
                ])
            restraints.append(lightpocket_restraint)

        
        hcdr3contact_restraint = ",".join([
                heavy_chain_letter, hcdr3_residx, antigen_chain_letter, antigen_residx, "contact",
                confidence, min_distance_angstrom, max_distance_angstrom, f"antibodycH-antigen"
            ])
        restraints.append(hcdr3contact_restraint)

    # Add restraint ids
    nr_restraints = len(restraints)
    restraints = [f"{restraints[j]},restraint_{j}" for j in range(nr_restraints)]

    # Write restraint output file
    restraint_out = "\n".join([chai_restraint_header] + restraints)
    with open(restraint_outfile, "w") as outfile:
        outfile.write(restraint_out)



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
### IMPORTS ###
from pathlib import Path
MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)

### FUNCTIONS ###
def abag_make_pocket_restraints(epitope_residues, restraint_outfile, antibody_chain_letters, confidence="1.0",
                                min_distance_angstrom="0.0", max_distance_angstrom="10.0"):
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

    connection_type = "pocket"
    chai_restraint_header = "chainA,res_idxA,chainB,res_idxB,connection_type,confidence,min_distance_angstrom,max_distance_angstrom,comment,restraint_id"

    # Collect antigen restraints
    restraints = []
    for j, epi in enumerate(epitope_residues):

        res, antigen_chain_letter, residx = epi
        antigen_residx = f"{res}{residx}"

        # antibody chain restraints 
        abchain1_restraint = ",".join([
            antibody_chain_letters[0], "", antigen_chain_letter, antigen_residx, connection_type,
            confidence, min_distance_angstrom, max_distance_angstrom, "antibodyc1-antigen"
        ])
        abchain2_restraint = ",".join([
            antibody_chain_letters[1], "", antigen_chain_letter, antigen_residx, connection_type,
            confidence, min_distance_angstrom, max_distance_angstrom, "antibodyc2-antigen"
        ])

        # Collect restraints
        restraints.append(abchain1_restraint)
        restraints.append(abchain2_restraint)

    # Add restraint ids
    nr_restraints = len(restraints)
    restraints = [f"{restraints[j]},restraint_{j}" for j in range(nr_restraints)]

    # Write restraint output file
    restraint_out = "\n".join([chai_restraint_header] + restraints)
    with open(restraint_outfile, "w") as outfile:
        outfile.write(restraint_out)
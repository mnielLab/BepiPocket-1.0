### IMPORTS ###

from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, Selection, PDBIO
import subprocess
from pathlib import Path
import pdb
import sys
MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)
AA3to1_DICT = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

### STATIC VARIABLES ###

ROOT_DIRECTORY = Path( Path(__file__).parent.resolve() )
MODELS_DIRECTORY = ROOT_DIRECTORY / "models"

### FUNCTIONS ###

def is_pdb_file(file_path):
    """
    File extension needs to be pdb/PDB, as well as readable by PDB parser.
    """

    # check file ending first
    if Path(file_path).suffix != '.pdb': return False
    # check if can be parsed as pdb, if no exception is raised, likely valid pdb file
    try:
        parser = PDBParser(QUIET=True)  
        structure = parser.get_structure('pdb_structure', file_path)
        return True  
    except Exception: return False
  
def is_cif_file(file_path):
    """
    File extension needs to be cif/CIF, as well as readable by CIF parser.
    """

    # check file ending first
    if Path(file_path).suffix != '.cif': return False
    # check if can be parsed as cif, if no exception is raised, likely valid cif file
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('cif_structure', file_path)
        return True

    except Exception: return False

def read_pdb_structure(pdb_file, pdb_id="foo", modelnr=0, return_all_models = False):
        """
        pdb_id: PDB acession, string
        pdb_file: path to 

        """
        #reading model 0 by default
    
        assert isinstance(modelnr, int), f"Model number needs to be a valid integer, it was {modelnr}"
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_file)

        #return all models
        if return_all_models:
            models = list()
            for m in structure: models.append(m)
            return models

        #return only desired model
        else: return structure[modelnr]

def read_cif_structure(cif_file, pdb_id="foo", modelnr=0, return_all_models=False):
    """
    Reads a CIF file and returns a specific model or all models in the structure.
    
    cif_file: path to the CIF file.
    pdb_id: PDB accession, string (default is 'foo').
    modelnr: The model number to return (default is 0).
    return_all_models: If True, returns all models in the structure.
    
    Returns the specified model or all models if return_all_models is True.
    """
    assert isinstance(modelnr, int), f"Model number needs to be a valid integer, it was {modelnr}"
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id, cif_file)

    if return_all_models:
        models = list(structure.get_models())
        return models
    else:
        return structure[modelnr]
    
def cif_to_pdb(cif_file, pdb_file, verbose=True):
    """
    Converts a CIF file to a PDB file using Biopython.

    Parameters:
        cif_file (str): Path to the input CIF file.
        pdb_file (str): Path to the output PDB file.

    Returns:
        None
    """
    try:
        # Parse the CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", cif_file)

        # Write the structure to a PDB file
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)

        if verbose: print(f"Successfully converted {cif_file} to {pdb_file}.")
    except Exception as e:
        print(f"An error occurred: {e}")

def write_biopdb_chain_residues_to_fasta(chains, pdb_acc_name, tgt_file=None):
    """
    Inputs: chains: List of Bio PDB chain class objects.
            epitope_or_paratope_res: List of Bio PDB residue class objects
    Outputs: AA_seqs. Amino acid sequences of chains.
    """
    #sometimes chains are passed as just one chain not in a list
    if not isinstance(chains, list): chains = [chains]
    
    fasta_file_content = str()
    AA_seqs = list()
    #remove hetero atoms
    for chain in chains: get_and_remove_heteroatoms(chain)

    for chain in chains:

        chain_id = chain.get_id()
        chain_residues = Selection.unfold_entities(chain, "R")
        
        AA_seq = str()

        for residue in chain_residues:
        
            try:
                aa = AA3to1_DICT[residue.get_resname()]
            #when residue is something nonstandard
            except KeyError:
                print(aa)
                print("Non-standard amino acid detected")
                aa = "X"
    
            AA_seq += aa

        #append to fasta_format
        fasta_file_content += f">{pdb_acc_name}_{chain_id}\n{AA_seq}\n"
        AA_seqs.append(AA_seq)
    if tgt_file != None:
        with open(tgt_file, "w") as outfile:
            outfile.write(fasta_file_content[:-1])

    return AA_seqs

def get_and_remove_heteroatoms(chain):
    """
   Heteroatoms in the form of water and other solvents need to be removed from the chain.
   Inputs: chain id

    """
    residues = Selection.unfold_entities(chain, "R")
    heteroatom_residue_ids = list()
    for residue in residues:
        residue_id = residue.get_full_id()

        #residue is a heteroatom
        if residue_id[3][0] != " ":
            heteroatom_residue_ids.append(residue_id[3])
    #remove all heteroatoms
    [chain.detach_child(ids) for ids in heteroatom_residue_ids]

def get_abag_interaction_data(ag_chains, ab_chains, return_bio_pdb_aas=False, atom_radius=4):
    
    epitope_data, paratope_data = [], []

    for ag_chain in ag_chains:
        for ab_chain in ab_chains:
            ab_epi_para_res = atom_neighbourhead_search_return_res(NeighborSearch( list(ag_chain.get_atoms()) ), list(ab_chain.get_atoms()), atom_radius=atom_radius)
            epitope_d, paratope_d = get_epitope_paratope_data(ab_epi_para_res, ag_chain, ab_chain, return_bio_pdb_aas=return_bio_pdb_aas)
            epitope_data.extend(epitope_d)
            paratope_data.extend(paratope_d)

    return epitope_data, paratope_data

def atom_neighbourhead_search_return_res(search_object, search_atoms, atom_radius=4):
    paired_interacting_residues = list()

    for search_atom in search_atoms:
        interact_res = search_object.search(search_atom.coord, radius = atom_radius, level="R")

        if interact_res:
            search_residue = search_atom.get_parent()
            if len(interact_res) == 1:
                paired_interacting_residues.append([interact_res[0], search_residue] )
            elif len(interact_res) > 1:
                for int_r in interact_res: paired_interacting_residues.append([int_r, search_residue])

    return paired_interacting_residues

def get_epitope_paratope_data(paired_residues, ag_chain, lc_or_hc, return_bio_pdb_aas = False):
    """
    
    """
    ag_id, ab_id = ag_chain.get_id(), lc_or_hc.get_id()
    antigen_residues = list(ag_chain.get_residues())
    antibody_residues = list(lc_or_hc.get_residues())


    epitope_data = []
    paratope_data = []

    ag_seq = write_pdb_res_to_seq(antigen_residues)
    lc_or_hc_seq = write_pdb_res_to_seq(antibody_residues)
    #create dict with indexes.
    ag_lookup = {antigen_residues[i]:i for i in range(len(antigen_residues))}
    ab_lookup = {antibody_residues[i]:i for i in range(len(antibody_residues))}

    for pair in paired_residues:
        epi_res = pair[0]
        para_res = pair[1]
        try:
            epi_res_idx = ag_lookup[epi_res] 
            para_res_idx = ab_lookup[para_res]
    
        except ValueError:
            sys.exit(f"Could not find paired residues in bio pdb residue list: {pair}")

        if not return_bio_pdb_aas:
            #get amino acid name'
            epi_res = AA3to1_DICT[epi_res.get_resname()]
            para_res = AA3to1_DICT[para_res.get_resname()]

        epitope_data.append((epi_res, ag_id, epi_res_idx) )
        paratope_data.append( (para_res, ab_id, para_res_idx) )
        
    return epitope_data, paratope_data

def write_pdb_res_to_seq(residues):
    """
    residues: Bio PDB residues
    """
    AA_seq = str()
    get_and_remove_heteroatoms(residues)

    for residue in residues:
        try:
            aa = AA3to1_DICT[residue.get_resname()]
        #when residue is something nonstandard
        except KeyError:
            print(aa)
            print("Non-standard amino acid detected")
            aa = "X"

        AA_seq += aa

    return AA_seq

def get_epitope_patch_residues(pred_structure_file, pred_epitope_residue, antigen_seqidxs=None, patch_angradius=10):

    if is_pdb_file(pred_structure_file): pred_structure = read_pdb_structure(pred_structure_file)
    elif is_cif_file(pred_structure_file): pred_structure = read_cif_structure(pred_structure_file)
    pred_structure_chains = list(pred_structure.get_chains())
    chain_ids = [c.get_id() for c in pred_structure_chains]
    nr_chains = len(pred_structure_chains)
    # define antigen + antibody chains
    if antigen_seqidxs is None: # assume antibody chains are always the last entries
        ag_chains = pred_structure_chains[:-2]
        ab_chains = pred_structure_chains[-2:]
    else:
        ag_chains = [pred_structure_chains[idx] for idx in antigen_seqidxs]
        ab_chains = [pred_structure_chains[i] for i in range(nr_chains) if i not in antigen_seqidxs]

    # get bio pdb search residue for predicted epitope residue
    pred_epitope_residue, pred_epitope_chainletter, pred_epitope_residue_idx  = pred_epitope_residue
    search_residue = list(pred_structure_chains[chain_ids.index(pred_epitope_chainletter)].get_residues())[pred_epitope_residue_idx]  

    # create atom list of antigen atoms
    search_atoms = []
    for ag_chain in ag_chains: search_atoms.extend( list(ag_chain.get_atoms()) )
    residue_pairs_within_atomradius = atom_neighbourhead_search_return_res(NeighborSearch( list( search_residue.get_atoms()) ), search_atoms, atom_radius=patch_angradius)

    search_residue_interacting_residues = []
    for residue_pair in residue_pairs_within_atomradius:
        _, residue = residue_pair
        chain_id = residue.parent.id  # Chain ID
        hetflag, res_id, icode = residue.id  # Residue ID (tuple: (hetflag, resseq, icode))
        res_name = residue.resname  # 3-letter amino acid name
        res_name = AA3to1_DICT[res_name]   
        search_residue_interacting_residues.append( (res_name, chain_id, res_id) )
    

    
    search_residue_interacting_residues = list( set(search_residue_interacting_residues) )


    return search_residue_interacting_residues 

def collect_epitope_contacts(pred_structure_files, antigen_seqidxs=None, epipara_aang_distance=4):

    pred_epitope_contacts, pred_structure_chains = [], []
    
    for pred_structure_file in pred_structure_files:

        # get predicte structure protein chains
        if is_pdb_file(pred_structure_file): pred_structure = read_pdb_structure(pred_structure_file)
        elif is_cif_file(pred_structure_file): pred_structure = read_cif_structure(pred_structure_file)
        pred_structure_chains = list(pred_structure.get_chains())
        nr_chains = len(pred_structure_chains)
        # define antigen + antibody chains
        if antigen_seqidxs is None: # assume antibody chains are always the last entries
            ag_chains = pred_structure_chains[:-2]
            ab_chains = pred_structure_chains[-2:]
        else:
            ag_chains = [pred_structure_chains[idx] for idx in antigen_seqidxs]
            ab_chains = [pred_structure_chains[i] for i in range(nr_chains) if i not in antigen_seqidxs]
        
        # get epitope and paratope residues
        epitope_data, paratope_data = get_abag_interaction_data(ag_chains, ab_chains, atom_radius=epipara_aang_distance)
        pred_epitope_contacts.append(epitope_data)

    return pred_epitope_contacts, pred_structure_chains
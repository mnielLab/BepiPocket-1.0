### IMPORTS ###
import sys
import pdb
from pathlib import Path
import subprocess
MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)


### FUNCTIONS ###

def read_fasta_entries(fasta_path):

    """
    Read fasta entries.

    Returns:
        list of tuples:
            [(header, sequence), ...]
    """

    entries = []
    header = None
    seq_lines = []

    with open(fasta_path, "r") as infile:

        for line in infile:

            line = line.strip()

            if line == "":
                continue

            if line.startswith(">"):

                if header is not None:
                    entries.append((header, "".join(seq_lines)))

                header = line[1:].strip()
                seq_lines = []

            else:
                seq_lines.append(line)

    if header is not None:
        entries.append((header, "".join(seq_lines)))

    return entries


def run_anarci(fasta_path, out_path, scheme="imgt"):

    """
    Run ANARCI on fasta file.
    Writes output to out_path.
    """

    cmd = [
        "ANARCI",
        "-i", str(fasta_path),
        "-o", str(out_path),
        "-s", scheme
    ]

    subprocess.run(cmd, check=True)


def parse_anarci_output(anarci_outfile):

    """
    Parse ANARCI output file.

    Returns:
        dict:
            {seq_id: [(residue, position, chain_type), ...]}
    """

    sequences = {}
    current_id = None

    with open(anarci_outfile, "r") as infile:

        for line in infile:

            line = line.strip()

            if line == "":
                continue

            if line == "//":
                current_id = None
                continue

            if line.startswith("#"):

                if line.startswith("# ANARCI"):
                    continue
                if line.startswith("# Domain"):
                    continue
                if line.startswith("# Most significant"):
                    continue
                if line.startswith("# Scheme"):
                    continue
                if line.startswith("#|"):
                    continue

                current_id = line[1:].strip()
                sequences[current_id] = []
                continue

            if current_id is None:
                continue

            parts = line.split()

            if len(parts) < 3:
                continue

            chain_type = parts[0]
            position = parts[1]
            residue = parts[2]

            if chain_type not in ["H", "K", "L"]:
                continue

            if residue == "-":
                continue

            sequences[current_id].append((residue, position, chain_type))

    return sequences


def identify_chain_type(sequence_entries):

    """
    Determine if sequence is heavy or light chain.

    Returns:
        "H" or "L"
    """

    chain_types = [entry[2] for entry in sequence_entries]

    if "H" in chain_types:
        return "H"
    else:
        return "L"


def extract_cdr_regions(sequence_entries, original_sequence, scheme="imgt"):

    """
    Extract CDR1, CDR2, CDR3.

    Returns:
        dict with CDR residue tuples

    Tuple format:
        (amino_acid, chain_type, residx)

    residx is the 0-based residue index in the original fasta sequence.
    """

    if scheme == "imgt":
        cdr_ranges = {
            "CDR1": (27, 38),
            "CDR2": (56, 65),
            "CDR3": (105, 117)
        }
    else:
        raise ValueError("Only IMGT scheme implemented")

    cdrs = {
        "CDR1": [],
        "CDR2": [],
        "CDR3": []
    }

    numbered_sequence = "".join([entry[0] for entry in sequence_entries])
    start_idx = original_sequence.find(numbered_sequence)

    if start_idx == -1:
        raise ValueError("ANARCI-numbered sequence not found in original sequence")

    for i in range(len(sequence_entries)):

        residue, position, chain_type = sequence_entries[i]
        pos_num = int("".join([c for c in position if c.isdigit()]))
        original_residx = start_idx + i

        if chain_type == "H":
            chain_type_out = "H"
        else:
            chain_type_out = "L"

        for cdr_name, (start, end) in cdr_ranges.items():
            if start <= pos_num <= end:
                cdrs[cdr_name].append((residue, chain_type_out, original_residx))

    return cdrs


def anarci_extract_cdrs(fasta_path, outdir, scheme="imgt"):

    """
    Full pipeline:
        - run ANARCI
        - parse output
        - identify chain type
        - extract CDRs
        - keep non-antibody sequences as None
    """


    if not outdir.is_dir(): outdir.mkdir(parents=True)
    anarci_outfile = outdir / "anarci_output.txt"

    run_anarci(fasta_path, anarci_outfile, scheme=scheme)

    fasta_entries = read_fasta_entries(fasta_path)
    parsed_sequences = parse_anarci_output(anarci_outfile)

    results = []

    for seq_id, seq in fasta_entries:

        if seq_id not in parsed_sequences or not parsed_sequences[seq_id]:
            results.append((seq_id, None))
            continue
        sequence_entries = parsed_sequences[seq_id]
        chain_type = identify_chain_type(sequence_entries)
        cdrs = extract_cdr_regions(sequence_entries, seq, scheme=scheme)

        results.append((seq_id, {
            "chain_type": chain_type,
            "CDR1": "".join([x[0] for x in cdrs["CDR1"]]),
            "CDR2": "".join([x[0] for x in cdrs["CDR2"]]),
            "CDR3": "".join([x[0] for x in cdrs["CDR3"]]),
            "CDR1_positions": cdrs["CDR1"],
            "CDR2_positions": cdrs["CDR2"],
            "CDR3_positions": cdrs["CDR3"]
        }))

    return results


def get_hcdr3_center_residue(fasta_path, outdir, structure_letters):
    # get cdr annotation (anarci)
    cdr_results = anarci_extract_cdrs(fasta_path, outdir)
    
    hcdr3_residues = None
    antigen_letters, light_chain_letter, heavy_chain_letter = [], [], []
    for i, cdr_result in enumerate(cdr_results):
        header, cdr_res = cdr_result
        if cdr_res is None:
            antigen_letters.append(structure_letters[i])
            continue
        
        chain_type = cdr_res["chain_type"]
        
        if chain_type == "L": light_chain_letter.append(structure_letters[i])
        elif chain_type == "H":
            heavy_chain_letter.append(structure_letters[i])
            hcdr3_residues = cdr_res["CDR3_positions"]

        else:
            antigen_letters.append(structure_letters[i])

    
    if light_chain_letter: light_chain_letter = light_chain_letter[0]
    else: light_chain_letter = None

    if heavy_chain_letter: heavy_chain_letter = heavy_chain_letter[0]
    else: heavy_chain_letter = None


    if heavy_chain_letter is None:
        raise ValueError("Heavy chain not found in fasta input")

    if hcdr3_residues is None or len(hcdr3_residues) == 0:
        raise ValueError("HCDR3 residues not found in fasta input")
    
    # get HCDR3 center residue
    center_idx = len(hcdr3_residues) // 2
    center_hcdr3_residue = hcdr3_residues[center_idx]
    aa, chain_type, residx = center_hcdr3_residue
    hcdr3_residue = (aa, heavy_chain_letter, residx + 1)

    return antigen_letters, light_chain_letter, heavy_chain_letter, hcdr3_residue
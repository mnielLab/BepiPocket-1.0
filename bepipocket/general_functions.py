### IMPORTS 

import sys
from pathlib import Path
import pickle
import shutil
import numpy as np

MODULE_DIR = str( Path( Path(__file__).parent.resolve() ) )
sys.path.append(MODULE_DIR)

### FUNCTIONS ###

def flatten(container):
    """
    Input: container: list or tuple.
    Output: A generator object
    """

    #go through each objects in container (list or tuple)
    for i in container:
        if isinstance(i, (list,tuple)):
            #if list or tuple, recursively yielding elements to the generator object
            for j in flatten(i): yield j
        #if not (just a single element), yield to generator object
        else: yield i

def split_list(input_list, num_sublists):
    """
    Split lists into a specified number of (num_sublists) list of lists
    input_list: list
    num_sublists: int
    """
    # Calculate the size of each sublist
    sublist_size = len(input_list) // num_sublists
    remainder = len(input_list) % num_sublists

    # Initialize variables
    start = 0
    end = 0
    sublists = []

    # Create sublists
    for i in range(num_sublists):
        start = end
        end += sublist_size + (1 if i < remainder else 0)
        sublists.append(input_list[start:end])

    return sublists

def load_pickle_file(infile_path):
    with open(infile_path, "rb") as infile:
        pickle_data = pickle.load(infile)
    return pickle_data

def copy_stuff(source_files, tgt, add_to_filename="", silence=False):
    """
    Move source files to tgt
    source_files: list of paths. Needs to be files
    tgt: files will be copied to  here.
    """
    tgt = Path(tgt)
    if not tgt.is_dir(): tgt.mkdir(parents=True)

    nr_source_files = len(source_files)
    for i in range(nr_source_files): 
        source_file = Path(source_files[i])
        destination = tgt / f"{source_file.stem}{add_to_filename}{''.join(source_file.suffixes)}"
     
        #if file doesls not already exist. copy it. 
        if not destination.is_file(): shutil.copy(source_file, destination)

        else:
            if not silence: print("File already existed. Skipping")

        if not silence: print(f"Parsed {i+1} / {nr_source_files} to ({tgt})")

def copy_stuff_v2(source_files, tgts, rm_source_files=False, add_to_filename=""):
    """
    Move source files to tgts
    source_files: list of paths. Needs to be files
    tgts: files will be copied to here.
    """

    nr_source_files = len(source_files)
    for i in range(nr_source_files): 
        source_file = Path(source_files[i])
        tgt = Path(tgts[i])
        destination = tgt / f"{source_file.stem}{add_to_filename}{''.join(source_file.suffixes)}"
     
        #if file doesls not already exist. copy it. 
        if not destination.is_file():
            if not tgt.is_dir(): tgt.mkdir(parents=True)
            print(f"Copied {source_file} to {destination}")
            shutil.copy(source_file, destination)

        else: print("File already existed. Skipping")

        if rm_source_files: source_file.unlink()

        print(f"Parsed {i+1} / {nr_source_files} to ({tgt})")

def copy_from_to_files(from_files, to_files, rm_after_copy=False):
    nr_files = len(from_files)
    for i in range(nr_files):

        from_file = Path(from_files[i])
        to_file = Path(to_files[i])

        #source file has to exist
        if not from_file.is_file():
            print(f"Source file {from_file} did not exist. Skipping")
            continue

        #copy if target file not already exists.
        if not to_file.is_file():
            if not to_file.parent.is_dir(): to_file.parent.mkdir(parents=True)
            shutil.copy(from_file, to_file)
            print(f"Copied {from_file.name} to {to_file} ({i+1} / {nr_files})")
            if rm_after_copy:
                from_file.unlink()
                print(f"Removed {from_file}") 

        else: print(f"Destination file with same name already exists: {to_file}. Skipping")


def _run_complete(p: Path) -> bool:
    if not p.is_dir():
        return False
    return (len(list(p.glob("*.cif"))) == 5) and (len(list(p.glob("*.npz"))) == 5)

def _wipe_dir(p: Path) -> None:
    for f in p.glob("*"):
        if f.is_file():
            f.unlink()
        else:
            import shutil
            shutil.rmtree(f)

def get_highest_confidence_structure(out_path):
    score_files = list(out_path.glob("scores*"))
    scores = [np.load(f) for f in score_files]
    scores = np.asarray([0.8 * s["iptm"][0] + 0.2 * s["ptm"][0] for s in scores])

    best_idx = np.argmax(scores)
    best_score = scores[best_idx]
    best_scorefile = score_files[best_idx]

    rank1_structure = best_scorefile.parent / f"{best_scorefile.stem.replace('scores', 'pred')}.cif"
    print(f"Structure with highest confidence {rank1_structure} with {best_score}")

    return rank1_structure
### IMPORTS AND STATIC STUFF ###
from bepipocket.normal_run import normal_run
from bepipocket.bepipocket_run import bepipocket_run
from bepipocket.randompocket_run import randompocket_run
from bepipocket.discopocket_run import discopocket_run
from pathlib import Path
import argparse 

### COMMAND LINE ARGUMENTS ###
parser = argparse.ArgumentParser("Create CHAI-1 models with different modelling mode")
parser.add_argument("-i", required=True, action="store", dest="fasta_file", type=Path, help="Input directory of fasta input files for Chai-1 to run.")
parser.add_argument("-ir", action="store", dest="restraint_file", default=None, type=Path, help="Input directory for restraint files for Chai-1 to run. Only required if running a restraint mode.")
parser.add_argument("-o", required=True, action="store", dest="out_dir", type=Path, help="Structure output directory")
#parser.add_argument("-pred", action="store", choices=["normal", "restraint", "bepipocket", "discopocket"], required=True, dest="pred", help="Chai-1 structure modelling mode to run.") #TODO
parser.add_argument("-pred", action="store", choices=["normal", "restraint", "bepipocket", "random"], required=True, dest="pred", help="Chai-1 structure modelling mode to run.")
parser.add_argument("-nr_runs", action="store", dest="nr_runs", type=int, default=6, help="Number of runs for Chai-1 structure modelling.")
parser.add_argument("-agscores", action="store", dest="agscores", default=None, type=Path, help="(Not required, will just compute on runtime, if not specified). Path to dict containing precomputed antigen sequence score (BepiPred-3.0 etc.): {FGKAJ...:array([0.4,,0.3,0.5,0.6,0.8...])..}.")
parser.add_argument("-hcdr3_mode", action="store_true", dest="hcdr3_mode", help="Create contact restraints between HCDR3 center reisidues and predicted BepiPocket or DiscoPocket residues.")
parser.add_argument("-hobohm_patchradius", action="store", dest="hobohm_patchradius", default=None, type=float, help="Spread epitope scoring with hobohm 1 patch algorithm.")
parser.add_argument("-msa_directory", action="store", dest="msa_directory", default=None, type=Path, help="Look for MSA .pqt files with sequence hash filenames mathcing query sequences in this directory.")

# set variables
args = parser.parse_args()
fasta_file = args.fasta_file
out_dir = args.out_dir
pred = args.pred
restraint_file = args.restraint_file
nr_runs= args.nr_runs

agscores = args.agscores
hcdr3_mode = args.hcdr3_mode
hobohm_patchradius = args.hobohm_patchradius
msa_directory = args.msa_directory

# chai-1 normal prediction mode 
if pred == "normal":
    normal_run(fasta_file, out_dir, overwrite_earlier_jobcontent=False, seeds=nr_runs, msa_directory=msa_directory)
# chai-1 restraint prediction mode (User defined restraints, as described in Chai-1 documentation)
elif pred == "restraint":
    normal_run(fasta_file, out_dir, restraint_file=restraint_file, overwrite_earlier_jobcontent=False, nr_runs=nr_runs, msa_directory=msa_directory)
# chai-1 bepipocket (use BepiPred-3.0 to guide antibody-epitope restraints)
elif pred == "bepipocket":
    bepipocket_run(fasta_file, out_dir, hcdr3_mode=hcdr3_mode, nr_runs=nr_runs, bp3_score_lookup=agscores, msa_directory=msa_directory, hobohm_patchradius=hobohm_patchradius)

elif pred == "discopocket":
    discopocket_run(fasta_file, out_dir, nr_runs=nr_runs, msa_directory=msa_directory)


# chai-1 randompocket (use random guidance for antibody-epitope restraint (baseline) )
elif pred == "random":
    randompocket_run(fasta_file, out_dir, nr_runs=nr_runs, random_score_lookup=agscores, msa_directory=msa_directory)




# chai-1 discopocket (use DiscoTope-3.0 to guide antibody-epitope restraints)
#TODO 
# elif pred == "discopocket":
#     abag_disco_score_key = out_dir.name.split("_discotopemap")[0]
#     discopocket_chairun(fasta_file, out_dir, patch_mode=patch_mode, nr_runs=seeds,
#                          discotope3_score_lookup = disco3scores, msa_directory=msa_directory, abag_disco_score_key=abag_disco_score_key)

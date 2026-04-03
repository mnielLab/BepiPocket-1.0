# BepiPocket

Current deep learning methods, such as AlphaFold and Chai, can often create antibody-antigen (AbAg) structures with high confidence and accuracy.
However, these methods often fail to predict the correct antibody binding site, placing the antibody incorrectly on the antigen and converging on repeatedly predicting the same redundant binding mode.
Varying seeds and diffusion samples does not guarantee that a diverse set of binding modes explored.
**BepiPocket** is a simple approach for integrating B-cell epitope prediction tools, **BepiPred-3.0** and **DiscoTope-3.0**, to guide antibody-epitope restraints during Chai-1 structure prediction.
This vastly increase diversity of explored antibody binding sites as well as allowing Chai-1 to find substantialy more accurate AbAg structures with high confidence.

### TODOS
* **TODO**: Upload MMseqs2 MSA Code. MSAs can be provide in the same format used as Chai-1 (.pqt). But are not automatically created on inference. Need to add this code.
* **TODO**: Merge in DiscoTope-3.0 package environment.

## License 
BepiPocket is a tool developed by the Health Tech section at Technical University of Denmark (DTU). The code and data can be used freely by academic groups for non-commercial purposes.
If you plan to use these tools for any for-profit application, you are required to obtain a license (contact Morten Nielsen, morni@dtu.dk).
* If you use BepiPred-3.0 (BepiPocket) to guide to create antibody-epitope pockets, please get a BepiPred-3.0 license.
* If you use DiscoTope-3.0 (DiscoPocket) to guide to create antibody-epitope pockets, please get a DiscoTope-3.0 license.
* The tool uses Chai-1 for structural inference of antibody-antigen complexes. Chai-1 is not open-source, and you are therefore also required to get license from Chai Discovery. 

## Graphical Abstract
![Screenshot](GraphicalAbstract.png)

## Installation 

##

### Create Conda Environment
This tool requires Python 3.10 It should work for later python versions as well.
```
$ conda create -n bepipocket python=3.10
$ conda activate bepipocket
$ conda install pip
```
### Install Pip Packages 
First, download requirements.txt file. Then,
```
$ pip install -r requirements.txt # install package dependencies
$ pip install git+https://github.com/mnielLab/BepiPocket-1.0.git
```
## Usage 
After installation, BepiPocket can be run with the script 'main.py'. 

### Inputs 

### Outputs

### Example
Using the example fasta file, this code snippet does 4 runs, each producing 5 structures.
The first run is without any antibody-epitope restraint, just using seed 0.
The second, third and and fourth runs are made using antibody-epitope restraints constructed with BepiPred-3.0. 
```
python main.py -i ./examples/2j88_ag_A_ab_L_H.fasta -o ./example_output/bepipredmap -pred bepipocket -nr_runs 4
```
## Cite
If you found this tool useful in your research, please cite:<br>
[Pocket Restraints Guided by B-Cell Epitope Prediction improves Chai-1 Antibody-Antigen Structure Modeling](https://doi.org/10.1101/2025.09.17.676770)

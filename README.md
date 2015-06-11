# DRCLUSTER
### Dimensionality Reduction Cluster

Clustering FASTA datasets using dimensionality reduction algorithms.

# Overview

### BLAST/jackhmmer all vs. all sequence comparison
- BLAST: Convert FASTA file to protein database and search using FASTA file as query
- jackhmmer: PSI-BLAST-like iterative search of a FASTA query against itself

Results are stored in a abstracted pairwise, symmetric, similarity matrix.


### SVD/t-SNE hybrid dimensionality reduction

Converts high dimensional similarity matrix into low-dimensional embedding while maintaining local manifold structure.

Embedding is stored as a two-dimensional numpy array.

### Plotting

Uses matplotlib to plot embedding with basic annotation tools.


## INSTALLATION

Ensure that the file structure for the src/ folder is as follows: (It should already be correct!)

- src
	- pyclust.py
	- lib/
		- annotation.py
		- grouper.py
		- initrun.py
		- jsonconv.py
		- mds_calc.py
		- plotter.py
		- readargs.py
		- results_parser.py
		- tsne/
		- tsne_calc.py

### Installation
- We recommend you set up a virtual environment and run the command ```pip install -r requirements.txt```
- Then run ```pip install tsne```. (The setup for tsne requires the numpy module)

- Alternatively, install all the packages in requirements.txt manually.


Requirements.txt contents:
```
Cython==0.22
matplotlib==1.4.3
scikit-learn==0.16.1
scipy==0.15.1
tsne==0.1.1
```

- Also, you must have a working installation of either BLAST or jackhmmer for all vs. all sequence alignment.

### Walkthrough for usage

Create a folder for output storage. The name of the folder *must* be the name of the the FASTA file. FASTA file name *must* be in the format **example.fas**
- ```mkdir proteins``` for FASTA file *proteins.fas*

Move FASTA file to this folder and create a *temp/* directory.
- ```mv proteins.fas proteins/```
- ```mkdir temp/```

Run all vs. all sequence comparison. BLAST output must be in table format under the name *results.out*. Jackhmmer output must be in table format under the name *tbl.hits*.
- ```$HMMERBIN/jackhmmer --tblout tbl.hits proteins.fas proteins.fas```

Navigate run pyclust.py from src/ directory with the following flags:
- ```-dir proteins``` name of data folder (e.g. proteins)
- ```-parse``` parse into similarity matrix
- ```-cluster``` cluster using SVD/t-SNE
- ```-plot``` plot using Matplotlib
- E.g. ```python /Path/to/src/pyclust.py -dir proteins -parse -cluster -plot```

# DRCLUSTER
### Dimensionality Reduction Cluster

Clustering FASTA datasets using dimensionality reduction algorithms.


## Overview

### BLAST+/jackhmmer all vs. all sequence comparison
- BLAST+: Convert FASTA file to protein database and search using FASTA file as query
- jackhmmer: PSI-BLAST-like iterative search of a FASTA query against itself

Results are stored in a abstracted pairwise, symmetric, similarity matrix.

### SVD/t-SNE hybrid dimensionality reduction

Converts high dimensional similarity matrix into low-dimensional embedding while maintaining local manifold structure.

Embedding is stored as a two-dimensional numpy array.

### Plotting

Uses matplotlib to plot embedding with basic annotation tools.


## Pip installation
- We recommend you set up a virtual environment and run the command `pip install -r requirements.txt`
- Then run `pip install tsne`. (The setup for tsne requires the numpy module)

- Alternatively, install all the packages in requirements.txt manually.


Required modules:
---------------------
Cython==0.22
matplotlib==1.4.3
scikit-learn==0.16.1
scipy==0.15.1
tsne==0.1.1
---------------------

- Also, you must have a working installation of either BLAST+ or jackhmmer for all vs. all sequence alignment.
- Only BLAST+ will work, legacy BLAST will not work!


## Walkthrough for minimal usage

Specify the path to the FASTA file using the -f flag. This is require
- Ex. `-f path/to/fasta.ff'

Select sequence alignment type using the -search flag. If no flag is specified, default is HMMER.
- Ex. `-search hmmer`

Either provide a path to HMMER or BLAST executables, or a path to an alignment file.
- Ex. `-align my/hmmer/results`
- Ex. `-bin my/hmmer/bin`
- Ensure that the results file or bin correspond to selected sequence alignment type

Matplotlib generated plot with `-plot` flag

Complete DRCluster command:
- Ex. `drcluster -search hmmer -bin path/hmmer/bin -plot`


## Using annotated FASTA files

To get coloring based on modelability and PFAM, include the `-a` flag, and format FASTA headers into the following: (note that a comparison to PDB database is required)

>*name*;*pfamID*;dom#;domlength;fulllength;domevalue;fullevalue;%len;modelability

- Ex. `>TR42;PF00240;1;63;120;3e-05;4e-4;88;mod`
















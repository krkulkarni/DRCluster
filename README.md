Python-TSNE
===========

Python library containing T-SNE algorithms.

Algorithms
----------

<<<<<<< HEAD
### Barnes-Hut-SNE
=======
### BLAST+/jackhmmer all vs. all sequence comparison
- BLAST+: Convert FASTA file to protein database and search using FASTA file as query
- jackhmmer: PSI-BLAST-like iterative search of a FASTA query against itself
>>>>>>> 5776ea233f4e5aabdef9c4aadacfa99983bee44b

A python ([cython](http://www.cython.org)) wrapper for [Barnes-Hut-SNE](http://homepage.tudelft.nl/19j49/t-SNE.html) aka fast-tsne.

<<<<<<< HEAD
I basicly took [osdf code](https://github.com/osdf/py_bh_tsne) and made it pip compilant.

The wrapper was successfully tested on OSX (10.6/10.7), Ubuntu (11.04) and Arch Linux.
=======
### SVD/t-SNE hybrid dimensionality reduction
>>>>>>> 5776ea233f4e5aabdef9c4aadacfa99983bee44b

Requirements
------------

* [numpy](numpy.scipy.org)>=1.7.1
* [scipy](http://www.scipy.org/)>=0.12.0
* [cython](cython.org)>=0.19.1
* [cblas](http://www.netlib.org/blas/) or [openblas](https://github.com/xianyi/OpenBLAS). Tested version is v0.2.5 and v0.2.6 (not necessary for OSX).

Installation
------------

You can install the package from [PyPI](https://pypi.python.org/pypi):

```
pip install tsne
```

<<<<<<< HEAD
Or directly from the Github repository:

```
pip install git+https://github.com/danielfrg/tsne.git
```
=======
## Pip installation
- We recommend you set up a virtual environment and run the command ```pip install -r requirements.txt```
- Then run ```pip install tsne```. (The setup for tsne requires the numpy module)
>>>>>>> 5776ea233f4e5aabdef9c4aadacfa99983bee44b

Usage
-----

Basic usage:

<<<<<<< HEAD
=======
Required modules:
>>>>>>> 5776ea233f4e5aabdef9c4aadacfa99983bee44b
```
from tsne import bh_sne
X_2d = bh_sne(X)
```

<<<<<<< HEAD
### Examples

* [Iris](http://nbviewer.ipython.org/urls/raw.github.com/danielfrg/py_tsne/master/examples/iris.ipynb)
* [MNIST](http://nbviewer.ipython.org/urls/raw.github.com/danielfrg/py_tsne/master/examples/mnist.ipynb)
* [word2vec on presidential speeches](https://github.com/prateekpg2455/U.S-Presidential-Speeches) via [@prateekpg2455](https://github.com/prateekpg2455)

More Information
----------------

See *Barnes-Hut-SNE* (2013), L.J.P. van der Maaten. It is available on [arxiv](http://arxiv.org/abs/1301.3342).
=======
- Also, you must have a working installation of either BLAST+ or jackhmmer for all vs. all sequence alignment.
- Only BLAST+ will work, legacy BLAST will not work!

## Walkthrough for minimal usage

Specify the path to the FASTA file using the -f flag. This is require
- ```-f path/to/fasta.ff```

Select sequence alignment type using the -search flag. If no flag is specified, default is HMMER.
- ```-search hmmer```

Either provide a path to HMMER or BLAST executables, or a path to an alignment file.
- ```-align my/hmmer/results```
- ```-bin my/hmmer/bin```
- Ensure that the results file or bin correspond to selected sequence alignment type

Matplotlib generated plot with ```-plot``` flag

Complete DRCluster command:
- ```drcluster -f path/to/fasta.ff -search hmmer -bin path/hmmer/bin -plot```

## Using annotated FASTA files

To get coloring based on modelability and PFAM, format FASTA headers into the following: (note that a comparison to PDB database is required)

> *>*name*;*pfamID;dom#;domlength;fulllength;domevalue;fullevalue;%len;modelability

- Ex. ```>TR42;PF00240;1;63;120;3e-05;4e-4;88;mod```














>>>>>>> 5776ea233f4e5aabdef9c4aadacfa99983bee44b


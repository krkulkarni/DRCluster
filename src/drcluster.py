#!/usr/bin/env python

import time
import sys
import os
import numpy as np
from lib import updated_readargs, initrun, results_parser, mds_calc, tsne_calc, grouper, plotter, runseqalign
from scipy import sparse
import tsne
#import json
#from lib import jsonconv

__author__ = "kulkarnik"

## The Algorithm class holds the three methods of dimensionality reduction.

## SVD - tSNE hybrid: svdsne()
## First performs singular value decomposition directly on sparse matrix
## to reduce the noise of dataset.
## Then, runs t-SNE with Barnes-Hut approximation on intermediate dataset
## to reduce to coordinate embedding

## Euclidean Multidimensional Scaling: mdsonly()
## Use with caution! It fails on large datasets
## Convert sparse matrix to dense matrix
## and perform linear dimensionality reduction
## treating similarity values as Euclidean distances

## t-SNE only: sneonly()
## Use with caution! This method fails under most conditions,
## due to the organization of the sparse similarity matrix
## Runs t-SNE algorithm directly on full matrix, without intermediate reduction,
## and without Barnes-Hut approximation

class Algorithm(object):

    def __init__(self,scipymat,pointslen,dim):
        self.scipymat = scipymat
        self.pointslen = pointslen
        self.dim = dim

    def svdsne(self,*args):
        print("Performing svdsne")
        tempred = min(self.pointslen/10,500)
        tempmatrix = mds_calc.svd(self.scipymat,tempred)
        matrix = tsne.bh_sne(tempmatrix)
        return matrix


    def mdsonly(self,*args):
        print("Performing mdsonly")
        if (self.pointslen > 2000):
            print "Too many proteins to perform MDS directly"
            sys.exit(2)
        matrix = mds_calc.metric_mds(self.scipymat.toarray(),self.dim)
        return matrix


    def sneonly(self,reinit,directory):
        print("Performing sneonly")
        inity = "{}/inity.npy".format(directory)
        if (reinit):
            try:
                os.remove(inity)
            except OSError:
                print "No initial inity file"
        if (self.pointslen > 2000):
            print "Too many proteins to perform t-SNE directly"
            sys.exit(2)
        matrix = tsne_calc.tsne(inity,False,scipymat.toarray(),
                                no_dims=self.dim,
                                initial_dims=self.pointslen)
        return matrix

class DRClusterRun(object):

    def __init__(self,args):
        # Store all arguments
        self.args = args
        self.base = self.args.fasta.split("/")[-1].split(".")[0]

        # Creates directory to hold parsed matrix and embedded coordinates (and other data)
        # Directory is created in the current working directory, with the base name of the fasta file
        if not (self.args.directory):
            if not os.path.exists("{}/".format(self.base)):
                os.makedirs(self.base)
            self.args.directory = "{}/".format(self.base)

        # Assign bin file based on alignment type
        if (self.args.search == 'blast'):
            self.exebin = self.args.blastpath
        elif (self.args.search == 'hmmer'):
            self.exebin = self.args.hmmerpath

        # Create dictionary of points from fastafile
        # See lib/initrun.py for more details
        self.points = initrun.read_fasta(args.fasta)


    def sequencealignment(self):
        # Run sequence alignment based on alignment type

        # Runobj is initialized with fastapath, path to executables,
        # and directory in which to store results
        runObj = runseqalign.Align(self.args.fasta,self.exebin,self.args.directory)

        if (self.args.search == 'hmmer'):
            ## Jackhmmer summary is stored in directory/base.jackhmmer_out
            jackhmmerout = "{}/{}.jackhmmer_out".format(self.args.directory,self.base)

            ## Jackhmmer all vs all output is stored in directory/base.jackhmmer_tbl
            output = runObj.runjackhmmer(5,1.0)
            self.args.alignfile = "{}/{}.jackhmmer_tbl".format(self.args.directory,self.base)
            with open(jackhmmerout, 'wb') as f:
                f.write(output[0][0])

        elif (self.args.search == 'blast'):
            ## BLAST all vs all output is stored in directory/base.blast_tbl
            output = runObj.runblast()
            self.args.alignfile = "{}/{}.blast_tbl".format(self.args.directory,self.base)

    def parseOutput(self):
        # Parse either the jackhmmer or BLAST output
        matrixpath = "{}/sparsedata.txt".format(self.args.directory)
        if not (self.args.preparsed):
            tabParser, tabHandle = initrun.open_file(self.args.alignfile)
            row,col,data = results_parser.next_line_original_format(self.args.value,
                                                                    tabParser,tabHandle,
                                                                    self.points,self.args.search)
            savemat = np.vstack((row,col,data))
            np.savetxt(matrixpath,savemat)
            scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(self.points),len(self.points)))
            return scipymat

        else:
            savemat = initrun.get_matrix(matrixpath)
            row = savemat[0]
            col = savemat[1]
            data = savemat[2]
            scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(self.points),len(self.points)))
            return scipymat

    def runClustering(self,scipymat):
        # Run the appropriate clustering algorithm
        # See lib/mds_calc.py for more details
        coordspath = "{}/{}_{}_coords.npy".format(self.args.directory,self.base,self.args.type)
        alg = Algorithm(scipymat,int(len(self.points)),int(self.args.dimension))

        options = {
            'svdsne': alg.svdsne,
            'mdsonly': alg.mdsonly,
            'sneonly': alg.sneonly
        }

        if not (self.args.preclustered):
            matrix = options[self.args.type](self.args.reinitialize,self.args.directory)
            np.savetxt(coordspath,matrix)

        else:
            matrix = np.loadtxt(coordspath)

        return matrix


    def createColors(self,matrix):
        # Initialize a colors array
        colors = np.zeros(shape=(len(self.points)))

        # Derive colors based on color flag
        if (self.args.color == 'pfam'):
            for name, point in self.points.iteritems():
                colors[point.index] = point.pfamnum
        elif (self.args.color == 'mod'):
            for name, point in self.points.iteritems():
                colors[point.index] = point.modcolor
        elif (self.args.color == 'group'):
            colors = grouper.dbscan(matrix)
            # colors = np.loadtxt(grouppath)
        else:
            print "Some error occurred"
            sys.exit(1)

        return colors

    def plotCoordinates(self,matrix,colors):
        sizespath = "{}/allsizes.txt".format(self.args.directory)
        if (self.args.directory.split('/')[-1] == "total_ub" or
                    self.args.directory.split('/')[-1] == "full_ub"):
            point_sizes = np.loadtxt(sizespath)
        else:
            point_sizes = 20
        if (int(self.args.dimension) == 2):
            plotter.pyplotter2d(matrix,colors,self.args.directory,self.points,point_sizes)

        elif (int(self.args.dimension) == 3):
            plotter.pyplotter3d(matrix,colors)

##RUN THE CODE
if (__name__ == '__main__'):
    t0 = time.time()
    print("Running script")

    args = updated_readargs.arg_parser()

    print("Initializing DRCluster")
    runclust = DRClusterRun(args)

    if not (args.alignfile):
        print("Do you want to run {} sequence alignment? y or n".format(args.search))
        if (raw_input() == 'y'):
            runclust.sequencealignment()
        else:
            print("Run stopped\n"
                  "Either run sequence alignment "
                  "or use the -align flag to specify a {} output file".format(args.search))
            sys.exit(0)

    print("Parsing output")
    scipymat = runclust.parseOutput()

    print("Running clustering algorithm")
    matrix = runclust.runClustering(scipymat)

    print ("Creating colors")
    colors = runclust.createColors(matrix)

    print("{} points in dataset".format(len(matrix)))
    print("Took {} seconds".format(time.time()-t0))

    print("Plotting with Matplotlib")
    runclust.plotCoordinates(matrix,colors)

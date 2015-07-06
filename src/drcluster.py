#!/usr/bin/env python -u

import time
import sys
import os
import csv
import numpy as np
from lib import updated_readargs, initrun, results_parser, grouper, plotter, runseqalign, algorithms
from scipy import sparse
#import json
#from lib import jsonconv

__author__ = "kulkarnik"
__version__ = "0.9.2"

class DRClusterRun(object):

    def __init__(self,args):
        # Store all arguments
        self.args = args
        self.base = self.args.fasta.split("/")[-1].split(".")[0]

        # Creates directory to hold parsed matrix and embedded coordinates (and other data)
        # Directory is created in the current working directory, with the base name of the fasta file
        try:
            if not (self.args.directory):
                if not os.path.exists("{}/".format(self.base)):
                    os.makedirs(self.base)
                self.args.directory = "{}/".format(self.base)

            elif not (os.path.exists(self.args.directory)):
                os.makedirs(self.args.directory)
        except OSError as err:
            print("Directory with FASTA name exists\n"
                  "Using directory with name {}_dir/".format(self.base))
            dirname = "{}_dir/".format(self.base)
            if not (os.path.exists(dirname)):
                os.makedirs(dirname)
            self.args.directory = dirname
        # Create dictionary of points from fastafile
        # See lib/initrun.py for more details
        self.points = initrun.read_fasta(self.args.fasta,self.args.annotated)


    def sequencealignment(self):
        # Run sequence alignment based on alignment type

        # Runobj is initialized with fastapath, path to executables,
        # and directory in which to store results
        runObj = runseqalign.Align(self.args.fasta,self.args.exebin,self.args.directory)

        if (self.args.search == 'hmmer'):
            ## Jackhmmer summary is stored in directory/base.jackhmmer_out
            jackhmmerout = "{}/{}.jackhmmer_out".format(self.args.directory,self.base)

            ## Jackhmmer all vs all output is stored in directory/base.jackhmmer_tbl
            try:
                output = runObj.runjackhmmer(5,args.evalue)
            except OSError:
                print("\nUnable to perform jackhmmer sequence alignment.\n"
                      "Check the commands above for correctness.\n"
                      "Did you pass the right executable directory with -bin?")
                sys.exit(1)
            self.args.alignfile = "{}/{}.jackhmmer_tbl".format(self.args.directory,self.base)
            with open(jackhmmerout, 'wb') as f:
                f.write(output[0][0])

        elif (self.args.search == 'blast'):
            ## BLAST all vs all output is stored in directory/base.blast_tbl
            try:
                output = runObj.runblast(args.evalue)
            except OSError:
                print("\nUnable to perform BLAST sequence alignment.\n"
                      "Check the commands above for correctness.\n"
                      "Did you pass the right executable directory with -bin?")
                sys.exit(1)
            self.args.alignfile = "{}/{}.blast_tbl".format(self.args.directory,self.base)

    def parseOutput(self):
        # Parse either the jackhmmer or BLAST output
        matrixpath = "{}/sparsedata.txt".format(self.args.directory)
        if not (self.args.preparsed):
            with open(self.args.alignfile) as f:
                totaloutputlen = sum(1 for _ in f)
            # tabParser, tabHandle = initrun.open_file(self.args.alignfile)
            row,col,data = results_parser.next_line_original_format(self.args.alignfile,
                                                                    self.args.value,
                                                                    self.points,self.args.search,
                                                                    totaloutputlen)
            savemat = np.vstack((row,col,data))
            print("Saving matrix data to {}".format(matrixpath))
            np.savetxt(matrixpath,savemat)
            print("Generating sparse matrix from data")
            scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(self.points),len(self.points)))
            return scipymat

        else:
            savemat = np.loadtxt(matrixpath)
            row = savemat[0]
            col = savemat[1]
            data = savemat[2]
            print("Generating sparse matrix from data")
            scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(self.points),len(self.points)))
            return scipymat

    def runClustering(self,scipymat):
        # Run the appropriate clustering algorithm
        # See lib/mds_calc.py for more details
        coordspath = "{}/{}_{}d_{}_coords.txt".format(self.args.directory,
                                                      self.base,self.args.dimension,
                                                      self.args.type)
        alg = algorithms.Algorithm(scipymat,self.points,int(self.args.dimension))

        if not (self.args.preclustered):
            if (self.args.type == "svdsne"):
                matrix = alg.svdsne(int(self.args.perplexity),float(self.args.theta))
            elif (self.args.type == "mdsonly"):
                matrix = alg.mdsonly()
            else:
                matrix = alg.sneonly(self.args.reinitialize,self.args.directory)

            fastacoords = "{}/{}_{}d_{}_fastacoords.txt".format(self.args.directory,
                                                                self.base,self.args.dimension,
                                                                self.args.type)
            with open(fastacoords, 'wb') as f:
                for _, point in self.points.iteritems():
                    f.write("{}\t".format(point.line[1:]))
                    for coord in matrix[point.index]:
                        f.write("{}\t".format(coord))
                    f.write("\n")

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
            grouppath = "{}/groups.txt".format(self.args.directory)
            colors = grouper.dbscan(matrix)
            grouper.savegroups(grouppath,self.points,colors)
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
    print("Running DRCluster (version {})".format(__version__))

    args = updated_readargs.arg_parser()

    print("Initializing DRCluster")
    runclust = DRClusterRun(args)

    print ("{} search type selected!".format(args.search))

    if (args.exebin):
        print("Executing {} sequence alignment".format(args.search))
        runclust.sequencealignment()

    print("Parsing output")
    scipymat = runclust.parseOutput()

    print("Running clustering algorithm")
    matrix = runclust.runClustering(scipymat)

    print ("Creating colors")
    colors = runclust.createColors(matrix)

    print("{} points in dataset".format(len(matrix)))
    print("Took {} seconds".format(time.time()-t0))

    if (args.plot):
        print("Plotting with Matplotlib")
        runclust.plotCoordinates(matrix,colors)

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


class DRClusterRun(object):

    def __init__(self,args):
        # Store all arguments
        self.args = args
        self.base = self.args.fasta.split("/")[-1].split(".")[0]

        # Create temp directory to hold parsed matrix and embedded coordinates (and other data)
        if not (self.args.directory):
            if not os.path.exists("{}/".format(self.base)):
                os.makedirs(self.base)
            self.args.directory = "{}/".format(self.base)

        # Assign executable based on alignment type
        if (self.args.search == 'blast'):
            self.executable = self.args.blastpath
        elif (self.args.search == 'hmmer'):
            self.executable = self.args.hmmerpath

        self.points = initrun.read_fasta(args.fasta)


    def sequencealignment(self):
        # Run sequence alignment based on alignment type
        runObj = runseqalign.Align(self.args.fasta,self.executable,self.args.directory)
        if (self.args.search == 'hmmer'):
            jackhmmerout = "{}/{}.jackhmmer_out".format(self.args.directory,self.base)
            output = runObj.runjackhmmer(5,1.0)
            self.args.alignfile = "{}/{}.jackhmmer_tbl".format(self.args.directory,self.base)
            with open(jackhmmerout, 'wb') as f:
                f.write(output[0][0])
        elif (self.args.search == 'blast'):
            output = runObj.runblast()
            self.args.alignfile = "{}/{}.blast_tbl".format(self.args.directory,self.base)

    def parseOutput(self):
        # Parse either the jackhmmer or BLAST output
        matrixpath = "{}/data.txt".format(self.args.directory)
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
        coordspath = "{}/{}_{}_coords.npy".format(self.args.directory,self.base,self.args.type)

        options = {
            'svdsne': self._svdsne,
            'mdsonly': self._mdsonly,
            'sneonly': self._sneonly
        }

        if not (self.args.preclustered):
            matrix = options[self.args.type](scipymat)
            np.save(coordspath,matrix)

        else:
            matrix = np.load(coordspath)

        return matrix

    def _svdsne(self,scipymat):
        tempred = min(int(len(self.points)/10),500)
        tempmatrix = mds_calc.svd(scipymat,tempred)
        matrix = tsne.bh_sne(tempmatrix)
        return matrix


    def _mdsonly(self,scipymat):
        if (len(self.points) > 2000):
            print "Too many proteins to perform MDS directly"
            sys.exit(2)
        matrix = mds_calc.metric_mds(scipymat.toarray(),int(self.args.dimension))
        return matrix


    def _sneonly(self,scipymat):
        inity = "{}/inity.npy".format(self.args.directory)
        if (self.args.reinitialize):
            try:
                os.remove(inity)
            except OSError:
                print "No initial inity file"
        if (len(self.points) > 2000):
            print "Too many proteins to perform t-SNE directly"
            sys.exit(2)
        matrix = tsne_calc.tsne(inity,False,scipymat.toarray(),
                                no_dims=int(self.args.dimension),
                                initial_dims=len(self.points))
        return matrix

    def createColors(self,matrix):
        colors = np.zeros(shape=(len(self.points)))

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
        if (self.args.directory.split('/')[-1] == "total_ub" or self.args.directory.split('/')[-1] == "full_ub"):
            point_sizes = np.loadtxt(sizespath)
        else:
            point_sizes = 20
        if (int(self.args.dimension) == 2):
            plotter.pyplotter2d(matrix,colors,self.args.directory,self.points,point_sizes)

        elif (int(self.args.dimension) == 3):
            plotter.pyplotter3d(matrix,colors)

##RUN THE CODE
if (__name__ == '__main__'):
    t0 = time.clock()
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
    print("Took {} seconds".format(time.clock()-t0))

    print("Plotting with Matplotlib")
    runclust.plotCoordinates(matrix,colors)


# def main():
#     ## Get the starting time to measure the time of the run
#     print "Running script..."
#     t0 = time.clock()
#
#     ## Obtain all the info from the arguments passed
#     args = updated_readargs.arg_parser()
#     print "Parsed arguments"
#
#     if not (args.directory):
#         if not os.path.exists("temp/"):
#             os.makedirs("temp")
#         args.directory = "temp/"
#     ## Define the paths of BLAST results file and names file
#     if (args.search == 'blast'):
#         resultsfile = args.alignfile
#         hmmerfile=''
#     elif (args.search == 'hmmer'):
#         hmmerfile = args.alignfile
#         resultsfile=''
#     else:
#         "Some error occurred"
#         raise EnvironmentError
#
#     base = args.fasta.split("/")[-1].split(".")[0]
#     matrixpath = "{}/data.txt".format(args.directory)
#     coordspath = "{}/{}_{}_coords.npy".format(args.directory,base,args.type)
#     jsonpath = "{}/file.json".format(args.directory)
#     inity = "{}/inity.npy".format(args.directory)
#     grouppath = "{}/groups.txt".format(args.directory)
#     sizepath = "{}/sizes.txt".format(args.directory)
#     reppath = "{}/reps.txt".format(args.directory)
#
#     if (args.reinitialize):
#         try:
#             os.remove(inity)
#         except OSError:
#             print "No initial inity file"
#
#     ## Obtain colors and names from names file
#     points = initrun.read_fasta(args.fasta)
#     print "Read FASTA file"
#
#     ## Check if coordinates have already been mapped
#     if (args.load):
#         print "Loading coordinates from path"
#         matrix = np.load(args.load)
#
#
#     elif (args.parse):
#
#         ## Obtain the handle to the results file
#         if (args.search == 'blast'):
#             tabParser, tabHandle = initrun.open_file(resultsfile)
#         elif (args.search == "hmmer"):
#             tabParser, tabHandle = initrun.open_file(hmmerfile)
#         else:
#             print("Some error occurred")
#             raise EnvironmentError
#
#         print "Opened BLAST results file"
#
#         # hdfmat, mathandle = initrun.create_matrix(args.value,points,matrixpath)
#         # print "Initialized matrix"
#
#         print "Parsing results"
#         row,col,data = results_parser.next_line_original_format(args.value,tabParser,tabHandle,points,args.search)
#         savemat = np.vstack((row,col,data))
#         np.savetxt(matrixpath,savemat)
#
#         scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(points),len(points)))
#
#     else:
#         print "Retrieving matrix"
#         savemat = initrun.get_matrix(matrixpath)
#         row = savemat[0]
#         col = savemat[1]
#         data = savemat[2]
#
#         scipymat = sparse.coo_matrix((data,(row,col)),shape=(len(points),len(points)))
#
#
#     ## Run the appropriate dimensionality reduction algorithm
#     ## -mdsonly = metric MDS with sklearn's manifold package
#     ## -svdsne = preprocess to "points/10" dimensions with MDS, then t-SNE for reduced matrix
#     ## -sneonly = t-SNE with van der Maatan algorithm
#
#     if (args.cluster):
#         if (args.type == "mdsonly"):
#             print "Performing MDS"
#             if (len(points) > 2000):
#                 print "Too many proteins to perform MDS directly"
#                 sys.exit(2)
#             matrix = mds_calc.metric_mds(scipymat.toarray(),int(args.dimension))
#
#         elif (args.type == "svdsne"):
#             print "Performing t-SNE with MDS preprocessing"
#
#             ## Partially reduce dimensionality of HDF5 matrix to 1/10th of original size or maximum of 400
#             tempred = min(int(len(points)/10),500)
#             print "Preprocessing the data using SVD..."
#             print "Reducing to", tempred, "dimensions"
#
#             tempmatrix = mds_calc.svd(scipymat,tempred)
#             matrix = tsne.bh_sne(tempmatrix)
#             # matrix = tsne_calc.tsne(inity,False,tempmatrix,no_dims=int(args.dimension),initial_dims=tempred)
#
#         elif (args.type == "sneonly"):
#             print "CAUTION: Performing t-SNE on full dissimilarity matrix"
#
#             if (len(points) > 2000):
#                 print "Too many proteins to perform t-SNE directly"
#                 sys.exit(2)
#             matrix = tsne_calc.tsne(inity,False,scipymat.toarray(),no_dims=int(args.dimension),initial_dims=len(points))
#         else:
#             print "Some error occurred! Please try again."
#             raise EnvironmentError
#
#         # save coordinates to file
#         np.save(coordspath,matrix)
#
#     else:
#         try:
#             print "Loading coordinates from temp: " + args.type
#             matrix = np.load(coordspath)
#         except IOError:
#             print "No coordinates in temp"
#
#     if (args.group):
#         labels = grouper.dbscan(matrix)
#         #groupids = grouping.findkmeans(matrix)
#         # print "Grouping points"
#         # groupids, reps, sizes = grouper.findgroups(matrix,args.group,points)
#         # np.savetxt(grouppath,groupids)
#         # np.savetxt(sizepath, sizes)
#         # grouper.savereps(reppath,reps,args.group)
#
#     # select correct color array
#     colors = np.zeros(shape=(len(points)))
#     if (args.color == 'pfam'):
#         for point in points:
#             colors[points[point].index] = points[point].pfamnum
#     elif (args.color == 'mod'):
#         for point in points:
#             colors[points[point].index] = points[point].modcolor
#     elif (args.color == 'group'):
#         colors = labels
#         # colors = np.loadtxt(grouppath)
#     else:
#         print "Some error occurred"
#         sys.exit(1)
#
#     print len(matrix), "points in dataset"
#     print "Took", time.clock()-t0, "seconds"
#
#     #with open(jsonpath, 'w') as jsonout:
#     #    json.dump(jsonconv.jsonmaker(colors,lines,matrix,args.format), jsonout, indent=2)
#
#     ## Plot the results with matplotlib's PyPlot
#     if (args.plot):
#         print "Plotting",len(matrix), "points"
#         ## load sizes manually!!!
#         if (args.directory.split('/')[-1] == "total_ub" or args.directory.split('/')[-1] == "full_ub"):
#             point_sizes = np.loadtxt(args.directory + "/allsizes.txt")
#         else:
#             point_sizes = 20
#         if (int(args.dimension) == 2):
#             plotter.pyplotter2d(matrix,colors,args.directory,points,point_sizes)
#
#         elif (int(args.dimension) == 3):
#             plotter.pyplotter3d(matrix,colors)





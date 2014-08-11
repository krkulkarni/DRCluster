<<<<<<< HEAD
import time, sys, os
import numpy as np
from Tkinter import *
from lib import readargs, initrun, results_parser, mds_calc, tsne_calc, plotter



__author__ = "kulkarnik"


##RUN THE CODE
if (__name__ == '__main__'):

    ## Get the starting time to measure the time of the run
    print "Running script..."
    t0 = time.clock()



    ## Obtain all the info from the arguments passed
    args = readargs.arg_parser()
    print "Parsed arguments"


    """
    If necessary, run BLAST locally (try to run on a faster machine instead)
    #initrun.run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff)
    """


    ## Define the paths of BLAST results file and names file
    resultsfile = args.directory + "/results.out"
    fastafile = args.directory + "/names_all"
    matrixpath = args.directory + "/temp/mds.hdf5"
    coordspath = args.directory + "/temp/coords.npy"
    inity = args.directory + "/temp/inity.npy"

    if (args.reinitialize):
        try:
            os.remove(inity)
        except OSError:
            pass

    ## Obtain colors and names from names file
    names,colors = initrun.read_fasta(fastafile,args.format)
    print "Read FASTA file"
    ## Obtain the handle to the results file
    tabParser, tabHandle = initrun.open_file(resultsfile)
    print "Opened file"

    if args.precoordinated == True:
        print "Loading coordinates"
        matrix = np.load(coordspath)
    else:
        ## Check if results are preparsed
        ##
        ## If no, create an HDF5 formatted matrix and initialize
        ## Then, parse BLAST results and populate the HDF5 matrix
        ##
        ## If yes, then obtain populated matrix from file

        if args.preparsed:
            print "Retrieving matrix"
            hdfmat = initrun.get_matrix(matrixpath)
        else:
            hdfmat = initrun.create_matrix(args.value,names,matrixpath)
            print "Initialized matrix"

            print "Parsing results"
            if (args.format == 'mod'):
                results_parser.next_line_modified_format(args.value,tabParser,tabHandle,names,hdfmat)
            elif (args.format == 'orig'):
                results_parser.next_line_original_format(args.value,tabParser,tabHandle,names,hdfmat)



    ## Check if coordinates have already been mapped



        ## Run the appropriate dimensionality reduction algorithm
        ## -mdsonly = metric MDS with sklearn's manifold package
        ## -snemds = preprocess to "points/10" dimensions with MDS, then t-SNE for reduced matrix
        ## -snepca = preprocess to "points/10" dimensions with PCA, then t-SNE for reduced matrix

        ## (still working on -n = nystrom MDS with pycogent's approximate_mds package)
        if (args.type == "mdsonly"):
            print "Performing MDS"
            matrix = mds_calc.metric_mds(hdfmat,int(args.dimension))

        elif (args.type == "snemds"):
            print "Performing t-SNE with MDS preprocessing"

            ## Partially reduce dimensionality of HDF5 matrix to 1/10th of original size or maximum of 400
            tempred = min(int(len(names)/10),400)
            print "Reducing to", tempred, "dimensions"

            tempmatrix = mds_calc.metric_mds(hdfmat,tempred)
            matrix = tsne_calc.tsne(inity,False,tempmatrix,no_dims=int(args.dimension),initial_dims=tempred)

        elif (args.type == "snepca"):
            print "Performing t-SNE with PCA preprocessing"

            ## Partially reduce dimensionality of HDF5 matrix to 1/10th of original size or maximum of 400
            tempred = min(int(len(names)/10),400)
            print "Reducing to", tempred, "dimensions"

            matrix = tsne_calc.tsne(inity,True,hdfmat,no_dims=int(args.dimension),initial_dims=tempred)

        elif (args.type == "sneonly"):
            print "CAUTION: Performing t-SNE on full dissimilarity matrix"

            if (len(names) > 2000):
                print "Too many proteins to perform t-SNE directly"
                sys.exit(2)
            matrix = tsne_calc.tsne(inity,False,hdfmat[...],no_dims=int(args.dimension),initial_dims=len(names))

        # model = TSNE(n_components=3,metric="precomputed")
        # matrix = model.fit_transform(hdfmat)
        #else: matrix = mds_calc.nystrom_frontend(len(names),math.sqrt(len(names)),2,mds_calc.getdist,hdfmat)

    # save coordinates to file
    np.save(coordspath,matrix)
    print "Took", time.clock()-t0, "seconds"


    ## Plot the results with matplotlib's PyPlot
    if (args.plot):
        print "Plotting",len(names), "points"

        if (int(args.dimension) == 2):
            plotter.pyplotter2d(matrix,colors,names)

        elif (int(args.dimension) == 3):
            plotter.pyplotter3d(matrix,colors)
=======
import math
import csv
import sys
import argparse
import h5py
import approximate_mds
import matplotlib.pyplot as plt
import time
from sklearn import manifold

##
##GLOBAL VARIABLES
##

__author__ = "kulkarnik"

resultsfile = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/jul_17_names/results.out"
fastafile = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/jul_17_names/names_new"

##
## Set up the argument parser using the argparse module
##

def arg_parser():
    paraParser = argparse.ArgumentParser(description='Clustering analysis of BLAST output using MDS algorithm')

    ## This mutually exclusive group chooses between the bit score and e-value
    ## as the value to generate the distance matrix
    valueToUse = paraParser.add_mutually_exclusive_group()
    valueToUse.add_argument('-b','--bitscore',help='Chooses the bitscore option', action="store_true", required = False)
    valueToUse.add_argument('-e','--evalue',help="Chooses the e-value option", action="store_true", required = False)

    ## This mutually exclusive group chooses between the 2D and 3D options to
    ## graph the protein clusters
    dimension = paraParser.add_mutually_exclusive_group()
    dimension.add_argument('-2d','--twod',help="Chooses the 2D plot function", action="store_true", required = False)
    dimension.add_argument('-3d','--threed', help="Chooses the 3D plot function", action="store_true", required = False)

    mdstype = paraParser.add_mutually_exclusive_group()
    mdstype.add_argument('-c','--classical',help="Chooses the classical MDS algorithm",action="store_true",required=False)
    mdstype.add_argument('-m','--metric',help="Chooses the non-metric MDS algorithm",action="store_true",required=False)


    args = paraParser.parse_args()

    ## These statements force user to choose one of either the -b or -e options,
    ## and one of either the -2d or -3d options
    if (args.bitscore == True):
        bitOrE = 'b'
    elif (args.evalue == True):
        bitOrE = 'e'
    else:
        print "Error, must pick -b or -e flag"
        sys.exit(2)

    if (args.twod== True):
        dim = 2
    elif (args.threed == True):
        dim = 3
    else:
        print "Error, must pick -2d or -3d flag"
        sys.exit(3)

    if (args.classical == True):
        mds = "c"
    elif (args.metric == True):
        mds = "m"
    else:
        print "Error, must pick -classical or -metric flag"
        sys.exit(3)

    ## return arguments for later use
    return bitOrE, dim, mds

##Read FASTA file names and make a list
def read_fasta(fastafile):
    names = []
    colors = []
    with open(fastafile) as f:
        for line in f:
            parts = line.split(";")
            names.append(parts[0].strip().split()[0])
            try:
                colors.append(int(parts[1]))
            except:
                pass
    return names,colors


##Create the distance matrix
##Initialize with 4 (the farthest possible value)
def create_matrix(flag):
    f = h5py.File("mds.h5","w")
    if (flag == 'b'):
        dset = f.create_dataset("hdfmat.h5",shape=(len(names),len(names)),fillvalue=4)
    else:
        dset = f.create_dataset("hdfmat.h5",shape=(len(names),len(names)),fillvalue=1)
    return dset


##
"""# TAB OPTION"""
##

##Creates handle for results.out file
##Parse tab delimited file to generate iterator
def open_file(filename):
    tabHandle = open(filename,"rb")
    tabParser = csv.reader(tabHandle, delimiter='\t')

    return tabParser, tabHandle

##Read each line in tab-delimited file and store important variables

##
## qSeqId   --> name of query
## qLen     --> length of query
## sSeqId   --> name of match
## sLen     --> length of match
## eValue   --> e-value of match
## bitScore --> bit score of match
##

def next_line(flag, parser, handle):
    ##if the bit flag is on, run addtoBitMatrix
    if (flag == 'b'):
        try:
            while (True):
                line = next(parser)
                qSeqId = line[0].split(";")[0]
                qLen = int(line[1].strip())
                sSeqId = line[2].split(";")[0]
                sLen = int(line[3].strip())
                bitScore = float(line[5].strip())
                add_to_bit_matrix(qSeqId,qLen,sSeqId,bitScore,sLen)

        except StopIteration:
            handle.close()

    ##otherwise run add to Ematrix
    else:
        try:
            while (True):
                line = next(parser)
                qSeqId = line[0].split(";")[0]
                qLen = int(line[1].strip())
                sSeqId = line[2].split(";")[0]
                eValue = float(line[4].strip())

                ##REMEMBER TO ADD FLAG OPTION FOR EITHER EVALUE OR BITSCORE MATRIX
                add_to_e_matrix(qSeqId,qLen,sSeqId,eValue)

        except StopIteration:
            handle.close()



##add scaled score to distance matrix
def add_to_bit_matrix(query,qLen,match,bits,sLen):

    ##look up query index and match index
    query_index = names.index(query)
    match_index = names.index(match)

    ##convert bit score into a scaled score
    bit_scaled_score = convert_bit_score(bits,qLen,sLen)

    #print query ,match, bit_scaled_score

    ##Only add scaled score to matrix if it is less than default and any other comparison
    if (hdfmat[query_index,match_index]>=4):
        hdfmat[query_index,match_index] = bit_scaled_score
        hdfmat[match_index,query_index] = bit_scaled_score


def convert_bit_score(bitscore,querylength,matchlength):
    divisor = min(querylength,matchlength)
    value = (math.log(.25/(bitscore/divisor))+2)
    return abs(value)

##WORK ON THE SCALED SCORE FOR E VALUES
def add_to_e_matrix(query,qLen,match,e):
    ##look up query index and match index
    query_index = names.index(query)
    match_index = names.index(match)

    ##convert bit score into a scaled score
    e_scaled_score = convert_e_score(e,qLen)

    #print query ,match, e_scaled_score

    ##Only add scaled score to matrix if it is less than default and any other comparison
    if (hdfmat[query_index,match_index]>=1):
        hdfmat[query_index,match_index] = e_scaled_score
        hdfmat[match_index,query_index] = e_scaled_score

def convert_e_score(evalue,querylength):
    if (evalue != 0):
        value = -math.log(evalue)
    else:
        value = 400

    if (value <= 1):
        value = 1

    value = (1/value)**0.3



    # evalue = evalue*1000
    # if (evalue < 0.000001):
    #     value = 0
    # elif (evalue > 10):
    #     value = 1
    # else:
    #     value = evalue/10
    return value

def getdist(i,j):
    return hdfmat[i,j]

##perform the MDS on HDF5 matrix
def cmds(matrix,dim):
    if (dim==2):
        ## MDS with 2 dimensions (default)
        (newmat,eig) = approximate_mds.cmds_tzeng(matrix,dim=2)
        ## call plotter with a 2D graph
        return newmat
    elif (dim==3):
        ## MDS with 3 dimensions
        (newmat,eig) = approximate_mds.cmds_tzeng(matrix,dim=3)
        ## call plotter with 3D graph
        return newmat

def nystrom_frontend(num_objects, num_seeds, dim, dist_func,permute_order=True):
    (seed_mat,restore_ids) = approximate_mds.build_seed_matrix(num_objects,num_seeds,dist_func,permute_order)
    mdscoords = approximate_mds.nystrom(seed_mat,dim)
    return mdscoords[restore_ids]

def metric_mds(mat):
    nmds = manifold.MDS(n_components=2, metric=True, max_iter=300,dissimilarity="precomputed", n_jobs=1,n_init=1,random_state=1)
    coords = nmds.fit(mat).embedding_
    return coords

def pyplotter2d(finalmat,mdstype,colors):
    if (mdstype == 'm'):
        x_points = finalmat[:,0]
        y_points = finalmat[:,1]
    elif (mdstype == 'c'):
        x_points = []
        y_points = []
        for pair in matrix:
            end = len(str(pair))
            x_points.append(str(pair)[2:end-2].strip().split()[0])
            y_points.append(str(pair)[2:end-2].strip().split()[1])


    fig, ax = plt.subplots()
    if (colors == []):
        colors = 'b'
    ax.scatter(x_points, y_points,c=colors)

    print len(names)
    print "Took", time.clock()-t0, "seconds"
    plt.show()

##RUN THE CODE
t0 = time.clock()
print "Running script..."
bitOrE, dim, mdstype = arg_parser()
print "Parsed arguments"

#run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff)
names,colors = read_fasta(fastafile)
print "Read FASTA file"

tabParser, tabHandle = open_file(resultsfile)
print "Opened file"

hdfmat = create_matrix(bitOrE)
print "Initialized matrix"

print "Parsing results"
next_line(bitOrE,tabParser,tabHandle)


print "Performing MDS"
if (mdstype == "m"):
    matrix = metric_mds(hdfmat)
elif (mdstype == "c"):
    matrix = cmds(hdfmat,dim)

#matrix = nystrom_frontend(len(names),math.sqrt(len(names)),2,getdist)

pyplotter2d(matrix,mdstype,colors)
>>>>>>> b49cf280b90621bc2f254f8cc63ae7e2bb2a0ebf






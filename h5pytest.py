import math
import csv
import sys
import argparse
import h5py
import approximate_mds
import matplotlib.pyplot as plt
import time


from Bio.Blast.Applications import NcbiblastpCommandline
import numpy
import rpy2.robjects as robj
from rpy2.robjects.packages import importr


plot3d = importr("scatterplot3d")
## for 3D plotting
MASS = importr("MASS")
## for non-metric MDS
rhdf5 = importr("rhdf5")
## for reading hdf matrix in R


##
##GLOBAL VARIABLES
##

__author__ = "kulkarnik"

resultsfile = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/newfas/results.out"
fastafile = "/Users/kulkarnik/Research/MDSCluster_2014/BLAST+/newfas/names_mainfam"

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
    mdstype.add_argument('-cmds','--classical',help="Chooses the classical MDS algorithm",action="store_true",required=False)
    mdstype.add_argument('-nmmds','--nonmetric',help="Chooses the non-metric MDS algorithm",action="store_true",required=False)


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
        mds = "cmds"
    elif (args.nonmetric == True):
        mds = "nmmds"
    else:
        print "Error, must pick -cmds or -nmmds flag"
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

    h5create = robj.r.h5createFile
    h5create("mds.h5")
    f = h5py.File("mds.h5","w")
    if (flag == 'b'):
        dset = f.create_dataset("hdfmat.h5",shape=(len(names),len(names)),fillvalue=1.0)

    else:
        dset = f.create_dataset("hdfmat.h5",shape=(len(names),len(names)),fillvalue=1.0)


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

    print query ,match, bit_scaled_score

    ##Only add scaled score to matrix if it is less than default and any other comparison
    if (bit_scaled_score < hdfmat[query_index,match_index] and bit_scaled_score < 3):
        hdfmat[query_index,match_index] = bit_scaled_score


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
    if (e_scaled_score < hdfmat[query_index,match_index] and e_scaled_score < 1):

        hdfmat[query_index,match_index] = e_scaled_score

def convert_e_score(evalue,querylength):
    if (evalue != 0):
        value = -math.log(evalue)
    else:
        value = 400

    if (value <= 1):
        value = 1

    value = (1/value)**0.4



    # evalue = evalue*1000
    # if (evalue < 0.000001):
    #     value = 0
    # elif (evalue > 10):
    #     value = 1
    # else:
    #     value = evalue/10
    return value

##convert numpy matrix to R matrix
def convert_to_r(mat):
    nr, nc = mat.shape
    mat_vec = robj.FloatVector(mat.transpose().reshape((mat.size)))
    r_mat = robj.r.matrix(mat_vec, nrow=nr, ncol=nc)
    return r_mat


##perform the MDS on HDF5 matrix
def mds(matrix,dim):

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


def point_plotter_2d(x,y):
    ##Define the R functions
    plot = robj.r.plot
    text = robj.r.text
    identify = robj.r.identify

    ## color array maker
    if (colors!=[]):
        groups = groupify(colors)
    else:
        groups = "black"
    ##Plot and label points
    plot(x,y, xlab='', ylab='',pch=16, col=groups)
    #identify(x,y,labels=names,cex=0.6,pos=4)
    #text(x, y, labels=names, cex=0.7, pos=4, col="black")

    ##Wait for user input to end
    raw_input()

def point_plotter_3d(x,y,z):
    ##Define the R functions
    plot = robj.r.scatterplot3d
    colors = robj.r.rainbow(201)
    col2rgb = robj.r.col2rgb
    t = robj.r.t
    text = robj.r.text
    identify = robj.r.identify


    ##Plot and label points
    groups = groupify(colors)
    plot(x,y,z, xlab='', ylab='',zlab='',pch=16)
    #identify(x,y,labels=names,cex=0.6,pos=4)
    #text(x, y, z, labels=names, cex=0.4, pos=4, col="black")

    ##Wait for user input to end
    raw_input()

def groupify(colors):

    c = robj.r.c

    rcolors = robj.IntVector(())
    for color in colors:
        rcolors = c(rcolors,color)

    return rcolors


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

print "Converting to numpy array"
## convert hdf5 array to numpy
numpyarr = hdfmat[:]

print "Performing MDS"
matrix = mds(hdfmat[...],dim)
print matrix

fig = plt.figure()
ax = fig.add_subplot(111)
x_points = []
y_points = []
for pair in matrix:
    end = len(str(pair))
    x_points.append(str(pair)[2:end-2].strip().split()[0])
    y_points.append(str(pair)[2:end-2].strip().split()[1])



p = ax.plot(x_points, y_points, 'ro')
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()

print "Took", time.clock()-t0, "seconds"
raw_input()

r_mat = convert_to_r(numpyarr)





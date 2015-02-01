__author__ = 'kulkarnik'
import argparse, os,sys

def arg_parser():
    paraParser = argparse.ArgumentParser(description='Clustering analysis of BLAST output using MDS algorithm')

    ## This argument accepts the name of the directory with all required files
    ##
    ## Inputted directory must contain:
    ## 1. Temp directory (this will store dissimilarity matrix [mds.hdf5] and calculated coordinates [coords.npy])
    ## 2. BLAST results file with the name "results.out"
    ## 3. Text file with namestrings of proteins with the name "names_all"
    ##      namestrings must be in the format: "name of protein","name of family" if you want coloring
    ##
    ## Clustering will create/read from:
    ## 1. temp/mds.hdf5
    ## 2. temp/coords.npy
    ## 3. temp/inity.npy (random generation of points for t-SNE clustering)
    curr = os.getcwd()

    paraParser.add_argument('-dir', '--directory',
                            help="Name of directory with all required info",
                            default=curr)

    ## This argument chooses between the bit score and e-value
    ## as the value to generate the distance matrix
    paraParser.add_argument('-val','--value',
                            help="Choose what value to use (bitscore or e-value)",
                            choices=['b','e'], default='e')

    ## This argument chooses between the 2D and 3D options to
    ## graph the protein clusters
    paraParser.add_argument('-dim','--dimension',
                            help="Choose how many dimensions to use",
                            choices=['2','3'], default='2')

    ## This argument chooses the type of clustering algorithm to use
    ## The choices are:
    ## snepca = run preprocessing of pairwise dissimilarity matrix with principal components analysis (PCA),
    ##          and run final clustering with t-SNE
    ## snemds = run preprocessing of pairwise dissimilarity matrix with multidimensional scaling (MDS),
    ##          and run final clustering with t-SNE
    ## mdsonly = run final clustering on pairwise dissimilarity matrix with MDS
    paraParser.add_argument('-type','--type',
                            help="Choose clustering algorithm",
                            choices=['snepca','snemds','mdsonly','sneonly'],default='snemds')

    ## This argument chooses the format of the BLAST results file
    ## For most results.out files, you should choose original format
    paraParser.add_argument('-fmat','--format',
                            help="In what format is BLAST output?",
                            choices=['orig','mod'], default='orig')

    ## Choose this argument if BLAST results are preparsed and HDF5 dissimilarity matrix has already been populated
    paraParser.add_argument('-pp', '--preparsed',
                            help="HDF-formatted distance matrix is already made",
                            action="store_true")

    paraParser.add_argument('-load', '--load',
                            help="Load coordinates from Numpy coordinate matrix")

    ## Choose this argument if clustering algorithm has already been run and coords.npy file is stored
    paraParser.add_argument('-pc', '--precoordinated',
                            help="MDS Calculation is already done",
                            action="store_true")

    ## Choose this argument to plot the coordinates in a PyPlot with matplotlib
    paraParser.add_argument('-plot','--plot',
                            help="Plot coordinates with matplotlib",
                            action="store_true")

    ## Choose this argument to create a new random initialization of points in t-SNE
    paraParser.add_argument('-reinit','--reinitialize',
                            help="Create new random initialization of points",
                            action="store_true")

    args = paraParser.parse_args()

    return args

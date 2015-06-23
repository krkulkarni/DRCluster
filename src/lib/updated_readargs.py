__author__ = 'kulkarnik'
import argparse, os,sys

def arg_parser():
    paraParser = argparse.ArgumentParser(description='Clustering analysis of BLAST output using MDS algorithm')

    paraParser.add_argument('-f', '--fasta',
                            help="Path to FASTA file",
                            required=True)

    paraParser.add_argument('-dir', '--directory',
                            help="Name of output directory",
                            default=None)

    ## This argument chooses between the bit score and e-value
    ## as the value to generate the distance matrix
    paraParser.add_argument('-val','--value',
                            help="Choose what value to use (bitscore or e-value)",
                            choices=['b','e'], default='e')

    paraParser.add_argument('-search', '--search',
                            help="Choose alignment program (BLAST or HMMER)",
                            choices=['blast', 'hmmer'],
                            default='hmmer')

    group = paraParser.add_mutually_exclusive_group(required=True)

    group.add_argument('-bin', '--exebin',
                            help="Path to binary executables, either BLAST or hmmer",
                            default=None)

    group.add_argument('-align', '--alignfile',
                            help="Jackhmmer or BLAST output file",
                            default=None)

    paraParser.add_argument('-e', '--evalue',
                            help="Evalue for sequence alignment",
                            type=float,default=1.0)

    ## This argument chooses between the 2D and 3D options to
    ## graph the protein clusters
    paraParser.add_argument('-dim','--dimension',
                            help="Choose how many dimensions to use",
                            choices=['2','3'], default='2')

    ## This argument chooses the type of clustering algorithm to use
    ## The choices are:
    ## svdsne = run preprocessing of (sparse) pairwise dissimilarity matrix with singular value decomposition (SVD),
    ##          and run final clustering with t-SNE
    ## mdsonly = run final clustering on pairwise dissimilarity matrix with MDS only
    ## sneonly = run final clustering on pairwise dissimilarity matrix with t-SNE only
    paraParser.add_argument('-type','--type',
                            help="Choose clustering algorithm",
                            choices=['svdsne','mdsonly','sneonly'],default='svdsne')

    paraParser.add_argument('-color','--color',
                            help="Choose coloring scheme: modelability, PFAM, or group",
                            choices=['pfam', 'mod', 'group'],default='mod')

    ## Choose this argument if BLAST results are preparsed and HDF5 dissimilarity matrix has already been populated
    paraParser.add_argument('-parsed', '--preparsed',
                            help="Results have already been parsed into a similarity matrix",
                            action="store_true")


    ## Choose this argument if clustering algorithm has already been run and coords.npy file is stored
    paraParser.add_argument('-clustered', '--preclustered',
                            help="Clustering algorithm has already been applied",
                            action="store_true")

    paraParser.add_argument('-group', '--group',
                            type=float,
                            help="Group into modeling families")

    ## Choose this argument to plot the coordinates in a PyPlot with matplotlib
    paraParser.add_argument('-plot','--plot',
                            help="Plot coordinates with matplotlib",
                            action="store_true")

    ## Choose this argument to create a new random initialization of points in t-SNE
    paraParser.add_argument('-reinit','--reinitialize',
                            help="Create new random initialization of points",
                            action="store_true")

    paraParser.add_argument('-a', '--annotated',
                            help="FASTA files are annotated with the PFAM and mod headers",
                            action="store_true")

    paraParser.add_argument('-perp', '--perplexity',
                            help="Set the number of neighbors to use in the t-SNE algorithm",
                            type=float,default=30.0)

    args = paraParser.parse_args()

    return args

__author__ = 'kulkarnik'
import argparse

def arg_parser():
    paraParser = argparse.ArgumentParser(description='Clustering analysis of BLAST output using MDS algorithm')

    ## This mutually exclusive group chooses between the bit score and e-value
    ## as the value to generate the distance matrix
    valueToUse = paraParser.add_argument('-val','--value',
                                         help="Choose what value to use (bitscore or e-value)",
                                         choices=['b','e'],required=True)

    ## This mutually exclusive group chooses between the 2D and 3D options to
    ## graph the protein clusters
    dimension = paraParser.add_argument('-dim','--dimension',
                                        help="Choose how many dimensions to use",
                                        choices=['2','3'], required = True)

    mdstype = paraParser.add_argument('-mds','--mdstype',
                                      help="Choose MDS algorithm",
                                      choices=['c','m'],required=True)

    fmat = paraParser.add_argument('-fmat','--format',
                                   help="In what format is BLAST output?",
                                   choices=['orig','mod'],required=True)

    pp = paraParser.add_argument('-pp', '--preparsed',
                                 help="HDF-formatted distance matrix is already made",
                                 action="store_true", required=False)

    precoords = paraParser.add_argument('-pc', '--precoordinated',
                                 help="MDS Calculation is already done",
                                 action="store_true", required=False)

    plot = paraParser.add_argument('-plot','--plot',
                                   help="Plot coordinates with matplotlib",
                                   action="store_true", required=False)
    dirname = paraParser.add_argument('-dir', '--directory',
                                      help="Name of directory with all required info",
                                      required=True)

    args = paraParser.parse_args()

    return args

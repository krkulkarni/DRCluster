__author__ = 'kulkarnik'
import argparse
import sys

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

    fmat = paraParser.add_mutually_exclusive_group()
    fmat.add_argument('-orig', '--original', help="Parse BLAST results in original format",action="store_true",required=False)
    fmat.add_argument('-mod', '--modified', help="Parse BLAST results in modified format",action="store_true",required=False)

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

    if (args.original == True):
        formatval = 'orig'
    elif (args.modified == True):
        formatval = 'mod'
    else:
        print "Error, must pick -original or -modified flag"
        sys.exit(4)
    ## return arguments for later use
    return bitOrE, dim, mds, formatval

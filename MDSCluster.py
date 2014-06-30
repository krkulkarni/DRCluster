from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import math
import csv
import numpy
import re
import sys, argparse
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
plot3d = importr("scatterplot3d")


##
##GLOBAL VARIABLES
##

__author__ = "kulkarnik"

resultsfile = "newfas/results.out"
fastafile = "newfas/names_mainfam"

queryname = "~/BLAST+/uniprot/uniprot.fas"
dbname = "~/BLAST+/uniprot/uniprot"
outfile = "~/BLAST+/uniprot/results.out"
fmt = '"6 qseqid qlen sseqid slen evalue bitscore"'
dbsize = 1000000
ecutoff = 1000000



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

    ## return arguments for later use
    return bitOrE, dim

##Set up the command for protein blast.
##Format 'blastp -query inputfile -db database -out tabfile -outfmt "6 qseqid qlen sseqid slen evalue bitscore" -dbsize 1000000'
##
##Results will be stored in an tab-delimited file in the chosen path
##Constant size of searchspace defined as 1,000,000 (TBD)
##
## CALL THIS FUNCTION IF RESULTS.OUT IS REMOVED
##

def run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff):
    blastP = NcbiblastpCommandline(query=queryname,
                                   db=dbname,
                                   out=outfile,
                                   outfmt=fmt,
                                   #dbsize=dbsize,
                                   #searchsp=10000
                                   )

    ##Run blastp locally and store results in results.out
    blastP()


##Read FASTA file names and make a list
def read_fasta(fastafile):
    names = []
    colors = []
    with open(fastafile) as f:
        for line in f:
            parts = line.split(";")
            names.append(parts[0].strip().split()[0])
            colors.append(int(parts[1]))
    return names,colors


##Create the distance matrix
##Initialize with 4 (the farthest possible value)
def create_matrix(flag):
    matrix = numpy.empty(shape=(len(names),len(names)))
    if (flag == 'b'):
        matrix.fill(500)
    else:
        matrix.fill(500)
    return matrix


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
    if (bit_scaled_score < matrix[query_index][match_index] and bit_scaled_score < 500):
        matrix[query_index][match_index] = bit_scaled_score


def convert_bit_score(bitscore,querylength,matchlength):
    divisor = min(querylength,matchlength)
    value = (math.log(.25/(bitscore/divisor))+2)
    return value

##WORK ON THE SCALED SCORE FOR E VALUES
def add_to_e_matrix(query,qLen,match,e):
    ##look up query index and match index
    query_index = names.index(query)
    match_index = names.index(match)

    ##convert bit score into a scaled score
    e_scaled_score = convert_e_score(e,qLen)

    print query ,match, e_scaled_score

    ##Only add scaled score to matrix if it is less than default and any other comparison
    if (e_scaled_score < matrix[query_index][match_index] and e_scaled_score < 500):
        matrix[query_index][match_index] = e_scaled_score

def convert_e_score(evalue,querylength):
    if (evalue!=0.0):
        value = math.log(evalue,2)
        if (value<-593.0):
            value = -593.0
    else:
        value = -593.0

    value = ((value+593.0)/100)**1.2
    return value

##convert numpy matrix to R matrix
def convert_to_r(mat):
    nr, nc = mat.shape
    mat_vec = robj.FloatVector(mat.transpose().reshape((mat.size)))
    r_mat = robj.r.matrix(mat_vec, nrow=nr, ncol=nc)
    return r_mat


##perform the MDS on R matrix
def mds(r_matrix,dim):
    ##Define the R functions
    cmd_scale = robj.r.cmdscale

    if (dim==2):
        ## MDS with 2 dimensions (default)
        points = cmd_scale(r_matrix)
        x = points.rx(True,1)
        y = points.rx(True,2)

        ## call plotter with a 2D graph
        point_plotter_2d(x,y)

    elif (dim==3):
        ## MDS with 3 dimensions
        points = cmd_scale(r_matrix,k=3)
        x = points.rx(True,1)
        y = points.rx(True,2)
        z = points.rx(True,3)

        ## call plotter with 3D graph
        point_plotter_3d(x,y,z)


def point_plotter_2d(x,y):
    ##Define the R functions
    plot = robj.r.plot
    text = robj.r.text
    identify = robj.r.identify

    ## color array maker
    groups = groupify(colors)

    ##Plot and label points
    plot(x,y, xlab='', ylab='',pch=16,col=groups)
    #identify(x,y,labels=names,cex=0.6,pos=4)
    #text(x, y, labels=names, cex=0.4, pos=4, col="black")

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
print "Running script..."
bitOrE, dim = arg_parser()
#run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff)
names,colors = read_fasta(fastafile)
tabParser, tabHandle = open_file(resultsfile)
matrix = create_matrix(bitOrE)
next_line(bitOrE,tabParser,tabHandle)
r_mat = convert_to_r(matrix)
mds(r_mat,dim)



##
"""XML OPTION"""
##

##
## REACTIVATE THESE LINES IF RESULTS.XML IS REMOVED
##
# def run_blast_xml():
#     blastXml = NcbiblastpCommandline(query="~/BLAST+/families/superfamilies.fas",
#                                    db="~/BLAST+/families/superfamilies",
#                                    out="results.xml",
#                                    outfmt=5,
#                                    dbsize=1000000
#                                    )
#
#
#     ##Run blastp locally and store results in results.xml
#     stdout, stderr = blastXml()
#

#
# ##Creates handle for results.xml file
# ##REMEMBER TO CLOSE HANDLE
# ##
# ##Parse XML file to generate iterator
# xmlhandle = open("results.xml")
# blast_records = NCBIXML.parse(xmlhandle)
#
#
# ##If enough space available, store results in python list.
# ##Try instead to iterate through, less space required
#
# #blast_records = list(blast_records)
#
#
# ##IMPORTANT VALUES:
#     ##record.query                        --> name of query protein
#     ##record.descriptions[i].e            --> E-value at i
#     ##record.query_letters                --> number of aa in query
#     ##record.descriptions[i].title        --> 'gn1|BL_ORD_ID|#' name of ith protein
#     ##
#     ##record.alignments[i].title          --> name as r.d[i].title
#     ##record.alignments[i].hsps[i].bits   --> bit score of ith protein
#     ##record.alignments[i].hsps[i].expect --> E-value at i
#
#
# ##Iterate through parsed records
# def nextRecord():
#
#     ## move to next record, and catch the StopIteration exception
#     try:
#         while (True):
#             record = next(blast_records)
#             addToMatrix(record)
#     except StopIteration:
#         handle.close()
#
#
# ##call function to add all values to matrix
# def addToMatrix(record):
#
#     ##store index of query
#     queryIndex = names.index(record.query)
#
#     ##for each alignment in the record,
#     ##remove BLAST tag and store e-value in matrix
#     for alignment in record.alignments:
#         for hsp in alignment.hsps:
#             match = re.sub(r'gnl\|BL_ORD_ID\|\d* ',r'',alignment.title)
#             matchindex = names.index(match)
#             matrix[queryIndex][matchindex] = hsp.expect
#             #print queryindex, ":::", matchindex, ":::", hsp.expect
#





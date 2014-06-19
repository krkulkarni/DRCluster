from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import math
import csv
import numpy
import re
import sys, argparse
import rpy2.robjects as robj

##
##GLOBAL VARIABLES
##

__author__ = "kulkarnik"
paraParser = argparse.ArgumentParser(description='This is the argument parser.')
group = paraParser.add_mutually_exclusive_group()
group.add_argument('-b','--bitscore',help='Chooses the bitscore option', action="store_true", required = False)
group.add_argument('-e','--evalue',help="Chooses the e-value option", action="store_true", required = False)
args = paraParser.parse_args()

if (args.bitscore == True):
    flag = 'b'
elif (args.evalue == True):
    flag = 'e'
else:
    print "Error, must pick -b or -e flag"
    sys.exit(2)

##Set up the command for protein blast.
##Format 'blastp -query inputfile -db database -out XMLfile -outfmt 5'
##
##Results will be stored in an XML file in the current working directory
##Constant size of database defined as 1,000,000 (TBD)

##
## REACTIVATE THESE LINES IF RESULTS.XML IS REMOVED
##
def runblastxml():
    blastxml = NcbiblastpCommandline(query="~/BLAST+/families/superfamilies.fas",
                                   db="~/BLAST+/families/superfamilies",
                                   out="results.xml",
                                   outfmt=5,
                                   dbsize=1000000
                                   )


    ##Run blastp locally and store results in results.xml
    stdout, stderr = blastp()


##Same command as above, stored in tab delimited file
##
## REACTIVATE THESE LINES IF RESULTS.OUT IS REMOVED
##

def runblasttab():
    blastp = NcbiblastpCommandline(query="~/BLAST+/families/superfamilies.fas",
                                   db="~/BLAST+/families/superfamilies",
                                   out="~/BLAST+/families/results.out",
                                   outfmt='"6 qseqid qlen sseqid slen evalue bitscore"',
                                   evalue = 1000000000,
                                   dbsize=1000000
                                   )


    ##Run blastp locally and store results in results.xml
    stdout, stderr = blastp()


##Read FASTA file names and make a list
def readfasta():
    names = []
    with open("families/names_superfamilies") as f:
        for line in f:
            names.append(line.strip())
    return names

##Create the distance matrix
##Initialize with 1 (the farthest possible value)

def createMatrix():
    matrix = numpy.empty(shape=(len(names),len(names)))
    matrix.fill(4)
    return matrix


##
"""# TAB OPTION"""
##

##Creates handle for results.out file
##Parse tab delimited file to generate iterator
tabhandle = open("families/results.out")
tabparser = csv.reader(tabhandle, delimiter='\t')


##Read each line in tab-delimited file and store important variables
def nextLine(flag):
    ##if the bit flag is on, run addtoBitMatrix
    if (flag == 'b'):
        try:
            while (True):
                line = next(tabparser)
                qseqid = line[0]
                qlen = int(line[1].strip())
                sseqid = line[2]
                slen = int(line[3].strip())
                evalue = float(line[4].strip())
                bitscore = float(line[5].strip())

                addtoBitMatrix(qseqid,qlen,sseqid,bitscore)

        except StopIteration:
            tabhandle.close()

    ##otherwise run add to Ematrix
    else:
        print "still working on e matrix"
        sys.exit(2)
        # try:
        #     while (True):
        #         line = next(tabparser)
        #         qseqid = line[0]
        #         qlen = int(line[1].strip())
        #         sseqid = line[2]
        #         slen = int(line[3].strip())
        #         evalue = float(line[4].strip())
        #         bitscore = float(line[5].strip())
        #
        #         ##REMEMBER TO ADD FLAG OPTION FOR EITHER EVALUE OR BITSCORE MATRIX
        #         addtoEMatrix(qseqid,qlen,sseqid,evalue)
        #
        # except StopIteration:
        #     tabhandle.close()



##add scaled score to distance matrix
def addtoBitMatrix(query,qlen,match,bits):
    ##look up query index and match index
    queryindex = names.index(query)
    matchindex = names.index(match)

    ##convert bit score into a scaled score
    bitscaledscore = (math.log(.25/(bits/qlen))+2)
    #print query ,match, bitscaledscore

    ##Only add scaled score to matrix if it is less than default and any other comparison
    if (bitscaledscore < matrix[queryindex][matchindex] and bitscaledscore < 4):
        matrix[queryindex][matchindex] = bitscaledscore



def addtoEMatrix(query,qlen,match,e):
    ##look up query index and match index
    queryindex = names.index(query)
    matchindex = names.index(match)

    ##convert bit score into a scaled score
    escaledscore = (math.log(.25/(e/qlen))+2)
    #print query ,match, escaledscore

    ##Only add scaled score to matrix if it is less than default and any other comparison
    if (escaledscore < matrix[queryindex][matchindex] and escaledscore < 4):
        matrix[queryindex][matchindex] = escaledscore


##convert numpy matrix to R matrix
def convertToR(mat):
    nr, nc = mat.shape
    matvec = robj.FloatVector(mat.transpose().reshape((mat.size)))
    rmat = robj.r.matrix(matvec, nrow=nr, ncol=nc)
    return rmat


##perform the MDS on R matrix
def mds(rmatrix):
    ##Define the R functions
    cmdscale = robj.r.cmdscale


    points = cmdscale(rmatrix)
    identify = robj.r.identify

    ##Obtain x-, y-coordinates from MDS
    x = points.rx(True,1)
    y = points.rx(True,2)

    return x,y

def pointPlotter(x,y):
    ##Define the R functions
    plot = robj.r.plot
    text = robj.r.text
    ##Plot and label points
    plot(x,y, xlab='', ylab='')
    #identify(x,y,labels=names,cex=0.6,pos=4)
    text(x, y, labels=names, cex=0.4, pos=4, col="black")

    ##Wait for user input to end
    raw_input()

##RUN THE CODE
runblasttab()
names = readfasta()
matrix = createMatrix()
nextLine(flag)
rmat = convertToR(matrix)
x,y = mds(rmat)
pointPlotter(x,y)


##
"""XML OPTION"""
##
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





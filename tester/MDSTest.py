from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import math
import numpy
import string
import re
import rpy2.robjects as robj

##Set up the command for protein blast.
##Format 'blastp -query inputfile -db database -out XMLfile -outfmt 5'
##
##Results will be stored in an XML file in the current working directory
##Constant size of database defined as 1,000,000 (TBD)
blastp = NcbiblastpCommandline(query="~/BLAST+/tester/oneprot.fas",
                               db="~/BLAST+/tester/oneprot",
                               out="results.xml",
                               outfmt=5,
                               dbsize=1000000
                               )


##Run blastp locally and store results in results.xml
stdout, stderr = blastp()

##Read FASTA file names and make a list
names = []
with open("name_oneprot") as f:
    for line in f:
        names.append(line.strip().strip('.'))


##Create the distance matrix
##Initialize with 1 (the farthest possible value)

matrix = numpy.ones(shape=(len(names),len(names)))

##Creates handle for results file
##REMEMBER TO CLOSE HANDLE
##
##Parse XML file to generate iterator
handle = open("results.xml")
blast_records = NCBIXML.parse(handle)


##If enough space available, store results in python list.
##Try instead to iterate through, less space required

#blast_records = list(blast_records)


##IMPORTANT VALUES:
    ##record.query                        --> name of query protein
    ##record.descriptions[i].e            --> E-value at i
    ##record.query_letters                --> number of aa in query
    ##record.descriptions[i].title        --> 'gn1|BL_ORD_ID|#' name of ith protein
    ##
    ##record.alignments[i].title          --> name as r.d[i].title
    ##record.alignments[i].hsps[i].bits   --> bit score of ith protein
    ##record.alignments[i].hsps[i].expect --> E-value at i


##Iterate through parsed records
def nextRecord():

    ## move to next record, and catch the StopIteration exception
    try:
        while (True):
            record = next(blast_records)
            addToMatrix(record)
    except StopIteration:
        handle.close()


##call function to add all values to matrix
def addToMatrix(record):

    ##store index of query
    query = record.query.encode("utf-8")
    queryindex = names.index(query)

    ##for each alignment in the record,
    ##remove BLAST tag and store e-value in matrix
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            match = re.sub(r'gnl\|BL_ORD_ID\|\d* ',r'',alignment.title).encode("utf-8")
            match = match.strip().strip(".")
            matchindex = names.index(match)
            matrix[queryindex][matchindex] = hsp.expect
            #print queryindex, ":::", matchindex, ":::", hsp.expect


##convert numpy matrix to R matrix
def convertToR(mat):
    nr, nc = mat.shape
    matvec = robj.FloatVector(mat.transpose().reshape((mat.size)))
    rmat = robj.r.matrix(matvec, nrow=nr, ncol=nc)
    return rmat

##perform the MDS on R matrix
def mds(rmatrix):
    cmdscale = robj.r.cmdscale
    plot = robj.r.plot
    loc = cmdscale(rmatrix)
    x = loc.rx(True,1)
    y = loc.rx(True,2)
    print x,y
    
nextRecord()
rmat = convertToR(matrix)
mds(rmat)





__author__ = 'kulkarnik'
import math
from scipy import sparse
import sys, time
import csv

##Read each line in tab-delimited file and store important variables

##
## qSeqId   --> name of query
## qLen     --> length of query
## sSeqId   --> name of match
## sLen     --> length of match
## eValue   --> e-value of match
## bitScore --> bit score of match
##

def next_line_original_format(alignfile, flag, points,search,totaloutputlen):
    ##if the bit flag is on, run addtoBitMatrix
    row = []
    col = []
    data = []
    print("Parsing line {} of {}".format("0",totaloutputlen))
    t0 = time.time()
    parser = csv.reader(open(alignfile,'rb'),delimiter='\t')

    for i, line in enumerate(parser,start=1):
        if (search == 'blast'):
            qSeqId = line[0].split(";")[0]
            sSeqId = line[1].split(";")[0]
            queryLen= line[3].split(";")[0]
            eValue = float(line[10].strip())
            bitscore = float(line[11].strip())
        elif (search == 'hmmer'):
            line = line[0].split()
            qSeqId = line[0].split(";")[0]
            sSeqId = line[2].split(";")[0]
            bitscore = float(line[5].strip())
            queryLen=100
            eValue = float(line[4].strip())
        else:
            print("Search flag error")
            sys.exit(1)
        ##REMEMBER TO ADD FLAG OPTION FOR EITHER EVALUE OR BITSCORE MATRIX
        if (flag == 'b'):
            row, col, data = add_to_bit_matrix(qSeqId,sSeqId,bitscore,queryLen,points,row,col,data)
        elif (flag == 'e'):
            row, col, data = add_to_e_matrix(qSeqId,sSeqId,eValue,points,row,col,data)
        if (i%50000 == 0):
            t1=time.time()-t0
            t0=time.time()
            print ("Parsing line {} of {} (50k lines in {} seconds)".format(i,totaloutputlen,t1))

    return row, col, data

##add scaled score to distance matrix
def add_to_bit_matrix(query,match,bitscore,querylen,points,row,col,data):

    ##look up query index and match index
    try:
        query_index = points[query].index
        match_index = points[match].index

    except KeyError as err:
        print err
        print("The protein name in alignment output life was not found!\n"
              "Are you sure that the FASTA file corresponds to the alignment file?")
        sys.exit(1)

    ##convert bit score into a scaled score
    bit_scaled_score = convert_bit_score(bitscore,querylen)

    ##Only add scaled score to matrix if current entry at position is default value

    if (query_index==match_index):
        return row,col,data

    row.append(query_index)
    col.append(match_index)
    data.append(bit_scaled_score)
    row.append(match_index)
    col.append(query_index)
    data.append(bit_scaled_score)

    return row,col,data


def convert_bit_score(bitscore,querylength):
    value = (math.log(.25/(bitscore/querylength))+2)
    return abs(value)

##WORK ON THE SCALED SCORE FOR E VALUES
def add_to_e_matrix(query,match,e,points,row,col,data):
    ##look up query index and match index
    try:
        query_index = points[query].index
        match_index = points[match].index

    except KeyError as err:
        print err
        print("The protein name in alignment output life was not found!\n"
              "Are you sure that the FASTA file corresponds to the alignment file?")
        sys.exit(1)
    ##convert bit score into a scaled score
    e_scaled_score = convert_e_score(e)

    ##Only add scaled score to matrix if current entry at position is default value
    if (query_index==match_index):
        return row,col,data

    row.append(query_index)
    col.append(match_index)
    data.append(e_scaled_score)
    row.append(match_index)
    col.append(query_index)
    data.append(e_scaled_score)

    return row,col,data

def convert_e_score(evalue):
    evalue = evalue/10
    try:
        if (evalue>1):
            return 0
        return -math.log(evalue)
        
    except ValueError:
        return 250

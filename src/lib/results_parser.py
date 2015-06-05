__author__ = 'kulkarnik'
import math
from scipy import sparse

##Read each line in tab-delimited file and store important variables

##
## qSeqId   --> name of query
## qLen     --> length of query
## sSeqId   --> name of match
## sLen     --> length of match
## eValue   --> e-value of match
## bitScore --> bit score of match
##
def next_line_original_format(flag, parser, handle,points,search):
    ##if the bit flag is on, run addtoBitMatrix
    row = []
    col = []
    data = []

    if (flag == 'b'):
        try:
            while (True):
                line = next(parser)
                qSeqId = line[0].split(";")[0]
                sSeqId = line[1].split(";")[0]
                qLen = int(line[7].strip())
                sLen = int(line[9].strip())
                bitScore = float(line[11].strip())
                add_to_bit_matrix(qSeqId,qLen,sSeqId,bitScore,sLen,points)

        except StopIteration:
            handle.close()

    ##otherwise run add to Ematrix
    else:
        try:
            if (search=='blast'):
                while (True):
                    line = next(parser)
                    qSeqId = line[0].split(";")[0]
                    sSeqId = line[1].split(";")[0]
                    eValue = float(line[10].strip())

                    ##REMEMBER TO ADD FLAG OPTION FOR EITHER EVALUE OR BITSCORE MATRIX
                    row, col, data = add_to_e_matrix(qSeqId,sSeqId,eValue,points,row,col,data)

            elif (search=='hmmer'):
                while (True):
                    line = next(parser)[0].split()
                    if (line[0].startswith("#")):
                        continue
                    qSeqId = line[0].split(";")[0]
                    sSeqId = line[2].split(";")[0]
                    eValue = float(line[4].strip())

                    ##REMEMBER TO ADD FLAG OPTION FOR EITHER EVALUE OR BITSCORE MATRIX
                    row, col, data = add_to_e_matrix(qSeqId,sSeqId,eValue,points,row,col,data)

        except StopIteration:
            handle.close()
            return row,col,data



##add scaled score to distance matrix
def add_to_bit_matrix(query,qLen,match,bits,sLen,points,hdfmat):

    ##look up query index and match index
    for i, item in enumerate(points):
        if item.name == query:
            query_index = i
        if item.name == match:
            match_index = i

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
def add_to_e_matrix(query,match,e,points,row,col,data):
    ##look up query index and match index

    query_index = points[query].index
    match_index = points[match].index

    # query_index = names.index(query)
    # match_index = names.index(match)

    ##convert bit score into a scaled score
    e_scaled_score = convert_e_score(e)

    #print query ,match, e_scaled_score

    ##Only add scaled score to matrix if current entry at position is default value

    # for r,c in zip(row,col):
    #     if not (r in row and c in col):
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

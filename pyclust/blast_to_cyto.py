__author__ = 'kulkarnik'
import csv

def open_file(filename):
    tabHandle = open(filename,"rb")
    tabParser = csv.reader(tabHandle, delimiter='\t')

    return tabParser, tabHandle

def next_line(parser, handle, ppfile, pptable):
    ##if the bit flag is on, run addtoBitMatrix

    try:
        while (True):
            line = next(parser)
            qSeqId = line[0].split(";")[0]
            sSeqId = line[1].split(",")[0]
            if (qSeqId == sSeqId):
                continue
            qLen = float(line[7].strip())
            sLen = float(line[9].strip())
            pIdent = float(line[2])
            pLenMatch = qLen/sLen

            if ((pIdent >= 30) and (pLenMatch >=0.97) and (pLenMatch <=1.03)):
                addToCyto(qSeqId, sSeqId, pIdent/100, pLenMatch, True, ppfile, pptable)
            else:
                addToCyto(qSeqId, sSeqId, pIdent/100, pLenMatch, False, ppfile, pptable)
    except StopIteration:
        print "Done!"
        handle.close()

def addToCyto(query, target, identity, lenmatch, interaction, ppfile, pptable):
    if (interaction and query!=target):
        ppfile.write(query+"\tpp\t"+target)
        ppfile.write("\n")
        pptable.write(query+","+target+","+str(identity)+","+str(lenmatch)+"\n")

if __name__ == "__main__":

    blastfile = "/Users/kulkarnik/Research/MDSCluster_2014/pyclust/ub50/results.out"
    ppfile = "/Users/kulkarnik/Research/MDSCluster_2014/pyclust/ub50/cyto.sif"
    pptable = "/Users/kulkarnik/Research/MDSCluster_2014/pyclust/ub50/cyto.csv"

    parser, handle = open_file(blastfile)
    with open(ppfile, 'w') as ppfile:
        with open(pptable, 'w') as pptable:
            next_line(parser,handle,ppfile,pptable)


## +"\t"+str(identity)+"\t"+str(lenmatch)
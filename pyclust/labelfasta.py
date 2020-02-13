__author__ = 'kulkarnik'

## USED TO LABEL FASTA WITH PDB SEQUENCE, MODELABLE, OR NON-MODELABLE

from Bio.Blast.Applications import NcbiblastpCommandline
import csv

def open_blastfile(filename):
    tabHandle = open(filename,"rb")
    tabParser = csv.reader(tabHandle, delimiter='\t')

    return tabParser, tabHandle


def checker(parser, handle):
    ##if the bit flag is on, run addtoBitMatrix
    status = ''
    ismodelable = False
    notmod = False
    try:
        while (True):
            line = next(parser)
            pident = float(line[2])
            length = float(line[3])
            qlen = float(line[4])
            plen = (length/qlen) * 100

            if (pident >= 99.0 and plen >= 99.0):
                status = 'pdb'
                break

            elif (pident >= 50.0 and plen >=75.0):
                status = 'mod'
                break

    except StopIteration:
        status = 'notmod'
        handle.close()

    return status


if __name__ == "__main__":

    directory = "/Users/kulkarnik/Research/MDSCluster_2014/01092015_ub_list/label/"
    fasta = directory+'ub.fas'
    tempfasta = directory+"temp.fas"
    database = directory+"pdb"
    output = directory+"results.out"
    moddedfasta = directory+"modded.fas"

    blastp_command = NcbiblastpCommandline(query=tempfasta, db=database, outfmt='"6 qseqid sseqid pident length qlen"', out=output)

    fastahandle = open(fasta, "rb")
    moddedhandle = open(moddedfasta, 'w')
    fastaparser = csv.reader(fastahandle,delimiter='\n')

    try:
        while (True):
            name = next(fastaparser)[0]
            seq = next(fastaparser)[0]
            with open(tempfasta, 'w') as t:
                t.write(name)
                t.write('\n')
                t.write(seq)

            stdout, stderr = blastp_command()
            parser, handle = open_blastfile(output)
            status = checker(parser,handle)

            name = name+";"+status
            moddedhandle.write(name)
            moddedhandle.write('\n')
            moddedhandle.write(seq)
            moddedhandle.write('\n')


    except StopIteration:
        print "\nfinished parsing"
    fastahandle.close()
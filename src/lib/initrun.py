__author__ = 'kulkarnik'
import sys
import h5py
import csv
import unittest
#from Bio.Blast.Applications import NcbiblastpCommandline

##Read FASTA file names and make a list

class AllPointInfo():

    def __init__(self, line, name, mod, modcolor, pfamnum, pfam, index):
        self.line = line
        self.name = name
        self.seq = None
        self.mod = mod
        self.modcolor = modcolor
        self.pfamnum = pfamnum
        self.pfam = pfam
        self.index = index

def read_fasta(fastafile):
    pfamdict = dict()
    colornum = 0
    points = {}
    with open(fastafile) as f:
        index = 0
        for line in f:
            if (line[0] == "#"):
                continue
            elif (line[0] == '>'):
                parts = line[1:].split(";")
                name = parts[0].strip().split()[0]
                modcolor = _convertmodtocolor(parts[8].strip())
                pfamnum, pfamdict, colornum = _checkpfam(parts[1].strip(),pfamdict,colornum)
                # create a new point
                #                      \full line   \name\ modelability  \mod color\ pfamnum\ pfam
                newpoint = AllPointInfo(line.strip(),name,parts[8].strip(),modcolor,pfamnum, parts[1].strip(),index)
                index = index+1

                sequence = f.next()
                newpoint.seq = sequence
                points[name] = newpoint

    return points

def _checkpfam(pfamname, pfamdict,colornum):
    if pfamname in pfamdict.keys():
        return pfamdict[pfamname], pfamdict, colornum
    else:
        colornum += 1
        pfamdict[pfamname] = colornum
        return colornum, pfamdict, colornum

def _convertmodtocolor(mod):
    if mod == 'pdb':
        return 1
    elif mod == 'mod':
        return 2
    elif mod == 'notmod':
        return 3

##Creates handle for results.out file
def open_file(filename):
    tabHandle = open(filename,"rb")
    tabParser = csv.reader(tabHandle, delimiter='\t')

    return tabParser, tabHandle


##Create the distance matrix
##Initialize with 4 (the farthest possible value)
def create_matrix(flag,points,matrixpath):
    f = h5py.File(matrixpath,"w")
    if (flag == 'b'):
        dset = f.create_dataset("dataset",shape=(len(points),len(points)),fillvalue=4)
    elif (flag == 'e'):
        dset = f.create_dataset("dataset",shape=(len(points),len(points)),fillvalue=0)
    return dset, f

def get_matrix(matrixpath):
    f = h5py.File(matrixpath,"r")
    dataset = f["dataset"]

    return dataset, f

# def run_blast_tab(queryname,dbname,outfile,fmt,dbsize,ecutoff):
#     blastP = NcbiblastpCommandline(query=queryname,
#                                    db=dbname,
#                                    out=outfile,
#                                    outfmt=fmt,
#                                    #dbsize=dbsize,
#                                    #searchsp=10000
#                                    )
#
#     ##Run blastp locally and store results in results.out
#     blastP()

##

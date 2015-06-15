__author__ = 'kulkarnik'
import sys
#import h5py
import numpy as np
import csv
from collections import OrderedDict
import unittest
#from Bio.Blast.Applications import NcbiblastpCommandline

##Read FASTA file names and make a list

class AllPointInfo():

    #def __init__(self, line, name, mod, modcolor, pfamnum, pfam, index):
    def __init__(self,line,name,index,
                 mod='mod',modcolor=2,
                 pfamnum=1,pfam='No fam'):
        self.line = line
        self.name = name
        self.seq = None
        self.mod = mod
        self.modcolor = modcolor
        self.pfamnum = pfamnum
        self.pfam = pfam
        self.index = index

def read_fasta(fastafile,annotated):
    pfamdict = dict()
    colornum = 0
    points = OrderedDict()
    with open(fastafile) as f:
        index = 0
        if (annotated):
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
                    newpoint = AllPointInfo(line.strip(),name, index,
                                            mod=parts[8].strip(),modcolor=modcolor,
                                            pfamnum=pfamnum, pfam=parts[1].strip())
                    index += 1

                    sequence = f.next()
                    newpoint.seq = sequence
                    points[name] = newpoint
        else:
            for line in f:
                if (line[0] == "#"):
                    continue
                elif (line[0] == ">"):
                    name = line[1:].strip().split(";")[0]
                    newpoint = AllPointInfo(line.strip(),name,index)
                    index += 1

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

# Retrieve matrix from path
def get_matrix(matrixpath):
    return np.loadtxt(matrixpath)

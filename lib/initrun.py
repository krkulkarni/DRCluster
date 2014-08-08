__author__ = 'kulkarnik'
import h5py
import csv
import unittest
#from Bio.Blast.Applications import NcbiblastpCommandline

##Read FASTA file names and make a list
def read_fasta(fastafile,fmat):
    names = []
    colors = []
    with open(fastafile) as f:
        for line in f:
            if (fmat == 'mod'):
                ## for modified format, namestrings are split at ";" character
                parts = line.split(";")
            elif (fmat == 'orig'):
                ## for original format, namestrings are split at "," character
                parts = line.split(",")
            names.append(parts[0].strip().split()[0])
            try:
                colors.append(int(parts[1]))
            except:
                pass
    return names,colors

##Creates handle for results.out file
def open_file(filename):
    tabHandle = open(filename,"rb")
    tabParser = csv.reader(tabHandle, delimiter='\t')

    return tabParser, tabHandle


##Create the distance matrix
##Initialize with 4 (the farthest possible value)
def create_matrix(flag,names,matrixpath):
    f = h5py.File(matrixpath,"w")
    if (flag == 'b'):
        dset = f.create_dataset("dataset",shape=(len(names),len(names)),fillvalue=4)
    else:
        dset = f.create_dataset("dataset",shape=(len(names),len(names)),fillvalue=1)
    return dset

def get_matrix(matrixpath):
    f = h5py.File(matrixpath,"r")
    dataset = f["dataset"]

    return dataset

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

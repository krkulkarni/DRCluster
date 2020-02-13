__author__ = 'kulkarnik'

import csv

def filterbyeval(fastafile,filterfile,cutoff):
    numofprots = 0

    with open(fastafile, 'rb') as f:
        with open (filterfile, 'w') as filter:
            fastaparser = csv.reader(f, delimiter='\n')

            try:
                while (True):
                    line = next(fastaparser)[0]
                    parts = line.split(";")
                    evalue = float(parts[5])
                    seq = next(fastaparser)[0]
                    if (evalue < cutoff):
                        filter.write(line)
                        filter.write('\n')
                        filter.write(seq)
                        filter.write('\n')
                        numofprots += 1

            except StopIteration:
                return numofprots

def filterbypfam(fastafile,filterfile,pfam):
    numofprots = 0

    with open(fastafile, 'rb') as f:
        with open (filterfile, 'w') as filter:
            fastaparser = csv.reader(f, delimiter='\n')

            try:
                while (True):
                    line = next(fastaparser)[0]
                    parts = line.split(";")
                    family = parts[1].split(".")[0]
                    seq = next(fastaparser)[0]
                    if (family == pfam):
                        filter.write(line)
                        filter.write('\n')
                        filter.write(seq)
                        filter.write('\n')
                        numofprots += 1
                return numofprots
            except StopIteration:
                return numofprots



if __name__ == "__main__":

    fasta = "/Users/kulkarnik/Research/MDSCluster_2014/01092015_ub_list/all_ub_labeled.fas"
    cutoff = 1.0e-37

    evaloutput = "/Users/kulkarnik/Research/MDSCluster_2014/01092015_ub_list/ub_" + str(cutoff) +".fas"


    # print "FASTA PATH: " + fasta
    # print "OUTPUT PATH: " + evaloutput
    #
    # numofprots = filterbyeval(fasta,evaloutput,cutoff)
    # print "NUMBER OF PROTEINS: " + str(numofprots)

    with open("/Users/kulkarnik/Research/MDSCluster_2014/01092015_ub_list/allfams", 'rb') as allfams:
        famparser = csv.reader(allfams, delimiter='\n')
        try:
            while (True):
                line = next(famparser)[0]
                print line
                pfamoutput = "/Users/kulkarnik/Research/MDSCluster_2014/01092015_ub_list/ub_" + line +".txt"
                numofprots = filterbypfam(fasta,pfamoutput,line)
                print "PFAM FILE: " + pfamoutput
                print "NUMBER OF PROTEINS: " + str(numofprots) + '\n'

        except StopIteration:
            print "Done"



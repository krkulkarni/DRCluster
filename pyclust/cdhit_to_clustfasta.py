__author__ = 'kulkarnik'
import csv

# accepts fasta and cluster files outputted from CD-HIT run (cluster and cluster.clstr)
def open_file(fastaname, clustername):
    fastaHandle = open(fastaname,"rb")
    clusterHandle = open(clustername,"rb")

    fastaParser = csv.reader(fastaHandle, delimiter='\n')
    clusterParser = csv.reader(clusterHandle, delimiter='\t')

    return fastaParser, clusterParser, fastaHandle, clusterHandle

def next_line(fastaparser, clusterparser, fastahandle, clusterhandle):
    # stores fas
    fastadict = {}
    repseq = {}
    try:
        while (True):
            line = next(fastaparser)[0]
            if (line[0] == ">"):
                name = line[1:]
                seq = next(fastaparser)[0]
                fastadict[name] = seq
    except StopIteration:
        print "Done with fasta file!"
        fastahandle.close()

    try:
        line = next(clusterparser)
        while (True):
            if (line[0][0] == ">"):
                pass
            else:
                id = line[1].split()
                name = id[1][1:len(id[1])-3]
                seq = fastadict.get(name)
                if (id[2] == "*"):
                    repseq[name] = seq
            line = next(clusterparser)

    except StopIteration:
        print "Done with cluster file!"
        clusterhandle.close()

    w = csv.writer(open(outputfasta, 'w'),delimiter='\n')
    for key, val in repseq.items():

        w.writerow([">"+key,val])

if __name__ == "__main__":

    fastafile = "/Users/kulkarnik/Research/MDSCluster_2014/cdhit/ubfam_125/cluster"
    clusterfile = "/Users/kulkarnik/Research/MDSCluster_2014/cdhit/ubfam_125/cluster.clstr"
    outputfasta = "/Users/kulkarnik/Research/MDSCluster_2014/cdhit/ubfam_125/repseq.fas"

    fastaparser, clusterparser, fastahandle, clusterhandle = open_file(fastafile, clusterfile)

    next_line(fastaparser,clusterparser,fastahandle,clusterhandle)


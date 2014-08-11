grep "^>" ~/BLAST+/families/superfamilies.fas | sed 's/^>//' > ~/BLAST+/families/names_superfamilies

python ~/BLAST+/MDSCluster.py -e -2d


# updated june 17, 2014 

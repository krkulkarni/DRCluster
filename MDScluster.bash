grep "^>" ~/BLAST+/families/superfamilies.fas | sed 's/^>//' > ~/BLAST+/families/names_superfamilies

python ~/BLAST+/MDSCluster.py

rm ~/BLAST+/families/names_superfamilies
rm results.xml

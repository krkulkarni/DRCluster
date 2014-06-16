from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import math

"""
Set up the command for protein blast.
Format 'blastp -query inputfile -db database -out XMLfile -outfmt 5'

Results will be stored in an XML file in the current working directory
Constant size of database defined as 1,000,000 (TBD)
"""

blastp = NcbiblastpCommandline(query="~/BLAST+/families/superfamilies.fas",
                               db="~/BLAST+/families/superfamilies",
                               out="results.xml",
                               outfmt=5,
                               dbsize=1000000
                               )
"""
Run blastp locally and store results in results.xml
"""

stdout, stderr = blastp()

"""
Creates handle for results file
REMEMBER TO CLOSE HANDLE

Parse XML file to generate iterator
"""

handle = open("results.xml")
blast_records = NCBIXML.parse(handle)

"""
If enough space available, store results in python list.
Try instead to iterate through, less space required
"""
## blast_records = list(blast_records)

"""
Iterate through parsed results
"""
record = next(blast_records)

"""
IMPORTANT VALUES:
record.descriptions[i].e            --> E-value at i
record.query_letters                --> number of aa in query
record.descriptions[i].title        --> 'gn1|BL_ORD_ID|#' name of ith protein

record.alignments[i].title          --> name as r.d[i].title
record.alignments[i].hsps[i].bits   --> bit score of ith protein

"""
for des in record.descriptions:
    value = (des.e/record.query_letters)/(1000000)
    value = math.log(value, 2)
    print des.title, value

print ""

for alignment in record.alignments:
    for hsp in alignment.hsps:
        value = (hsp.bits/record.query_letters)/2.08
        print alignment.title, value

handle.close()




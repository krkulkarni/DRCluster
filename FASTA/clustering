import csv
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import math

blastp = NcbiblastpCommandline(query="../families/superfamilies.fas", db="../families/superfamilies", out="results.xml", outfmt=5, dbsize=1000000)
stdout, stderr = blastp()

handle = open("results.xml")
blast_records = NCBIXML.parse(handle)

blast_records = list(blast_records)

record = blast_records[0]

for des in record.descriptions:
    value = (des.e/record.query_letters)/(1000000)
    value = math.log(value, 2)
    print (des.title, value)

print ""

for al in record.alignments:
    for hsp in al.hsps:
        print (al.title, hsp.score)






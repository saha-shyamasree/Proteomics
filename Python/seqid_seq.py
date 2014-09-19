from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re


handle = open(sys.argv[1], "rU")
records = list(SeqIO.parse(handle, "fasta"))
#print("hello")
for record in records:
    s_id=re.sub(r'^(?:sp|tr)\|(.*?)\|(.*)?',r'\1',record.id)
    print(s_id+"\t"+str(len(record.seq)))
handle.close()

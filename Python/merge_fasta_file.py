from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

for i in range(1,len(sys.argv)):
    handle = open(sys.argv[i], "rU")
    for record in SeqIO.parse(handle, "fasta"):
        print(">"+record.description)
        print(record.seq)
    handle.close()

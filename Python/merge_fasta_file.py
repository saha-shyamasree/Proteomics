from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

wrt=open(sys.argv[len(sys.argv)-1],"w")
for i in range(1,len(sys.argv)-1):
    handle = open(sys.argv[i], "rU")
    for record in SeqIO.parse(handle, "fasta"):
        wrt.write(">"+record.description+"\n")
        wrt.write(str(record.seq)+"\n")
    handle.close()
wrt.close()

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re

if len(sys.argv)!=3:
    raise Exception("Need one fasta and one id file name as argument")

with open(sys.argv[2],'rU') as f:
    overlapping_ids = list(f)
    overlapping_ids[:] = [re.sub(r'^(\S+)\s(.*)?',r'\1',line) for line in overlapping_ids]
    overlapping_ids[:] = [line.strip() for line in overlapping_ids]
"""
for line in overlapping_ids:
    #line=line.strip();
    print(line)
"""
handle = open(sys.argv[1], "rU")
records = list(SeqIO.parse(handle, "fasta"))
print("Total=")
print(len(records))
count1=0
count2=0
for record in records: # SeqIO.parse(handle, "fasta")
    #print(record.id)
    if record.id in overlapping_ids:
        count1=count1+1
        #print(">"+record.description)
        #print(record.seq)
        #continue
        #print("found:"+record.description)
    else:
        count2=count2+1
        #continue
        #print(">"+record.description)
        #print(record.seq)
    #print(record.id)
handle.close()
print("Count1="),
print(count1)
print("Count2="),
print(count2)
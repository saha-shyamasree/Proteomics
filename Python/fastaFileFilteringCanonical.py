from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re
import argparse
import csv

parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("-f","--fastaFile", help="fasta file name")
args = parser.parse_args()

#fastahandle = open(sys.argv[1], "rU")
#headerhandle = open(sys.argv[2], "rU")

fastahandle = open(args.fastaFile, "rU")
records = list(SeqIO.parse(fastahandle, "fasta"))
#print("Total=")
#print(len(records))
count1=0
count2=0
for record in records: # SeqIO.parse(handle, "fasta")
    if record.id.startswith("sp") and "-" not in record.id:
        count1=count1+1
        #print("id:"+record.id)
        print(">"+record.description)
        print(record.seq)
        #continue
        #print("found:"+record.description)
fastahandle.close()

#print(overlapping_ids)
print("Count:"+str(count1))
#C:\Users\shyama\Dropbox\results\Bat_human_hendra\Human\trinity_only.fasta
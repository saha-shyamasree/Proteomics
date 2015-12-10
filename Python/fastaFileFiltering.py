from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re
import argparse
import csv

parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("-f","--fastaFile", help="fasta file name")
parser.add_argument("-l","--chosenIds", help="Fasta header list file name, headers should be newline separated")
parser.add_argument("-o","--outFile", help="Output fasta file name")
args = parser.parse_args()

#fastahandle = open(sys.argv[1], "rU")
#headerhandle = open(sys.argv[2], "rU")

fastahandle = open(args.fastaFile, "rU")
headerhandle = open(args.chosenIds, "rU")
outputHandle = open(args.outFile, "w")

overlapping_ids = list(headerhandle)
overlapping_ids[:] = [line.strip() for line in overlapping_ids]
overlapping_ids[:] = [line.replace('"','') for line in overlapping_ids]
#print("oids:"+str(len(overlapping_ids)))
#print(overlapping_ids[0])
records = list(SeqIO.parse(fastahandle, "fasta"))
#print("Total=")
#print(len(records))
count1=0
count2=0
for record in records: # SeqIO.parse(handle, "fasta")
    #print(record.description)
    #print("TESTTEST")
    if record.description in overlapping_ids:
        overlapping_ids.remove(record.description)
        count1=count1+1
        outputHandle.write(">"+record.description+"\n")
        outputHandle.write(record.seq+"\n")
        #continue
        #print("found:"+record.description)
fastahandle.close()
headerhandle.close()
#print(overlapping_ids)
#print("Count:"+str(count1))
#C:\Users\shyama\Dropbox\results\Bat_human_hendra\Human\trinity_only.fasta
    
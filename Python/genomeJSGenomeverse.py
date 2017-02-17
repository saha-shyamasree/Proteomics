### This code reads genome fasta and creates javascript file required for Genoverse.

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re
import os
import argparse


def filterFasta(fastaFile, outFile, species):
	fastahandle = open(fastaFile, "rU")
	outputHandle = open(outFile, "w")
	records = list(SeqIO.parse(fastahandle, "fasta"))
	outputHandle.write("var "+species+" = {\n")
	for record in records:
		outputHandle.write("\t\""+record.id+"\": {\n")
		outputHandle.write("\t\t\"size\": "+str(len(record.seq))+",\n")
		outputHandle.write("\t\t\"bands\": []\n")
		outputHandle.write("\t},\n")
	outputHandle.write("};\n")
	fastahandle.close() 
	outputHandle.close()


parser = argparse.ArgumentParser(description='This script read genome.fasta and creates a javascript file to create karyotype of Genoverse')
parser.add_argument("-f", "--fasta", nargs=1, required=True, help="full path of the genome fasta file", metavar="PATH")
parser.add_argument("-s", "--species", nargs=1, required=True, help="Species name", metavar="PATH")

args = parser.parse_args()

print(args.fasta[0])

outFile = os.path.dirname(args.fasta[0])+"/"+args.species[0]+".js"
print(outFile)
filterFasta(args.fasta[0],outFile,args.species[0])

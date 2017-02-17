###this code reads annotation files from PITDB post processing and standard proteome, uniprot fasta file
###Take only the TGEs classified as isoform, novel, known with variations and take their uniprot map. These
###uniprot protein sequences are added to the PIT generated fasta file and spectra are searched against this
###fasta file. I need three parameters, annotation file, PIT fasta, standard fasta ()
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import re
import argparse
import pandas as pd

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def customFasta(orfsFasta, refsFasta, iso, outFile):
    orfsIds=iso['ORF Id'].tolist()
    ##merging the source, protein id and protein name because thats how the fasta id is constructed
    refsIds=(iso['Source']+"|"+iso['Protein ID']+"|"+iso['Protein Name']).tolist()
    #print(refsIds)
    orfHandle = open(orfsFasta, "rU")
    refsHandle = open(refsFasta, "rU")
    outHandle = open(outFile, "w")
    records = list(SeqIO.parse(orfHandle, "fasta"))
    count=0
    for record in records: # SeqIO.parse(handle, "fasta")
        if record.description in orfsIds:
            #print("found1")
            count=count+1
            outHandle.write(">"+record.description+"\n")
            outHandle.write(str(record.seq)+"\n")
    recordRef = list(SeqIO.parse(refsHandle, "fasta"))
    print("count="+str(count))
    count=0
    print(len(set(refsIds)))
    for record in recordRef: # SeqIO.parse(handle, "fasta")
        #print(record.id)
        #print(record.description)
        if record.id in refsIds:
            #print("found2")
            count=count+1
            outHandle.write(">"+record.description+"\n")
            outHandle.write(str(record.seq)+"\n")
    print("count="+str(count))
    orfHandle.close()
    refsHandle.close()
    outHandle.close()
    
parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("--ORFs", required=True, help="ORFs fasta file name")
parser.add_argument("--uniprots", required=True, help="Uniprot/reference fasta file name")
parser.add_argument("--annotation", required=True, help="annotation csv file")
parser.add_argument("--outFile", required=True, help="fasta out file name")

args = parser.parse_args()
annotation=readFile(args.annotation,',')
print(annotation.shape)
#isoforms=annotation[(annotation['Class']!="known") & (annotation['Class']!="novel")]
##isoforms including known variations, novel
#isoforms=annotation[annotation['Class']!="known"]
#All identified ORFs
isoforms=annotation
print(isoforms.shape)
customFasta(args.ORFs, args.uniprots, isoforms, args.outFile)
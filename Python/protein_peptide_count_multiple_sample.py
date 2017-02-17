## this code reads protein and peptide file and count total number of protein and peptides for a dataset, i.e. multiple samples.
## This code assume that the protein and peptide files of all the samples are in the same folder. Corresponding fatsa files are in another folder.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse
import glob
import pandas as pd
import re

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

##seq1 and seq2 are list of sequences.
def uniqueSequence(seq1, seq2):
    ##this function get two sets of fasta sequences and return unique list.
    seq=seq1 + seq2
    #set function is keeping only unique entries
    return list(set(seq))

def readFastaFile(fastaFile):
    ##This function reads a fasta file using SeqIO, convert entries into record object and put all the sequences together.
    fastaHandle = open(fastaFile, "rU")
    records = list(SeqIO.parse(fastaHandle, "fasta"))
    #print("Records:"+str(len(records)))
    seqList=[str(r.seq) for r in records]
    #print("seqList:"+str(len(seqList)))
    return seqList

 
parser = argparse.ArgumentParser(description='This python code counts unique number of proteins and peptides in a dataset')
parser.add_argument("-d", "--dir", nargs=1, required=True, help="full path of protein and peptide files folder", metavar="PATH")
parser.add_argument("-f", "--fasta", nargs=1, required=True, help="full path of fasta files folder", metavar="PATH")

args = parser.parse_args()
#print(args)
#onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
onlyfiles=glob.glob(args.dir[0]+"/*+fdr+th+grouping_filtered.csv")
#print(args.dir[0])
print(onlyfiles)
peptides=[]
proteins=[]
for f in onlyfiles:
    print("F:"+f)
    fBase=re.sub(re.escape("+fdr+th+grouping_filtered.csv"),"",f)
    #print("fBase"+fBase)
    sample=os.path.basename(fBase).split(".",maxsplit=1)[0]
    print("sample:"+str(sample))
    pepDF=readFile(f,',')
    peptides=peptides+pepDF['Sequence'].tolist()
    sProt=readFastaFile(args.fasta[0]+sample+".assemblies.fasta.transdecoder.pep.identified.fasta")
    proteins=uniqueSequence(proteins , sProt)

unqPeptides=list(set(peptides))
print("Proteins/TGEs:"+str(len(list(set(proteins)))))
print("Peptides/Sequences:"+str(len(unqPeptides)))
protFile=open(args.fasta[0]+"TGEs.tsv","w")
protFile.write("Sequence\n")
protFile.write("\n".join(proteins))
protFile.close()

pepFile=open(args.dir[0]+"uniqPeptides.tsv","w")
pepFile.write("Peptides\n")
pepFile.write("\n".join(unqPeptides))
pepFile.close()
    
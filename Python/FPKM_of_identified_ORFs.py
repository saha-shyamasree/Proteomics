## This code reads the RSEM stat file and filter that to keep only identified ones. One transcript
## may have several identified ORFs.
import os
import argparse
import pandas as pd
import glob
import re

def readFile(filename, sep, headerFlag):
	fileDFObj=None
	if os.path.getsize(filename)>0:
		if headerFlag==0:
		    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
		elif headerFlag==1:
		    fileDFObj = pd.read_table(filename, sep=sep, header=None, keep_default_na=False, na_values=[''])
		else:
		    print("Unrecognized Header Flag")
	return fileDFObj

def filter(fpkmFile, annotFile):
    ##this function read the RSEM produced isoform stat file and the annotation file to filter
    ##the RSEM file based on the annotation file entries.
    fpkm=readFile(fpkmFile, '\t', 0)
    annotation=readFile(annotFile,',',0)
    annotIds=annotation['ORF Id'].str.extract("([^\|]+)")
    print("Total Identified:"+str(annotation.shape[0]))
    print("Total Annot Ids:"+str(annotIds.size))
    identified=fpkm[fpkm['transcript_id'].isin(annotIds)]
    print("Total identified Ids:"+str(identified.shape[0]))    
    return identified
    
parser = argparse.ArgumentParser(description='This python code find failed transcripts for both the gmap and blat aligner. These are possible virus transcripts. This code then extracts the trinity trascript and save in a file for further processing.')
parser.add_argument("-f", "--fpkm", nargs=1, required=True, help="full path of the RSEM isoform stat file", metavar="PATH")
parser.add_argument("-a", "--annot", nargs=1, required=True, help="full path of the annotation file", metavar="PATH")
parser.add_argument("-o", "--out", nargs=1, required=True, help="full path of the annotation file", metavar="PATH")

args = parser.parse_args()
print(args)
identified=filter(args.fpkm[0], args.annot[0])
print(args.out[0])
identified.to_csv(args.out[0],sep="\t", index=False)
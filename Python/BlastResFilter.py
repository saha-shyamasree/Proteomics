##This code reads BLAST result in csv format (output of contigStat.py) and a list of ids (identified ORFs).
##Queries from the id list is written to a file.

import pandas as pd
import argparse

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def pepThresholding(prots, pepTh, protRevStr, protContStr):
    print("1. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protRevStr)]
    print("2. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protContStr)]
    print("3. prots dim:"+str(prots.shape))
    prots=prots[prots['distinct peptide sequences']>pepTh]
    print("4. prots dim:"+str(prots.shape))
    return prots

def filterBlast(blastFile, queryList, outFile):
    blast=readFile(blastFile,",")
    blast['query_name'].isin(queryList)
    identIdx=blast['query_name'].isin(queryList)
    idx=identIdx[identIdx==True].index.tolist()
    blastIdentified= blast.iloc[idx,:]
    blastIdentified.to_csv(outFile,index=False)
    
parser = argparse.ArgumentParser(description='This python code filters blast results based on identified ORFs.')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="BLAST result for all the ORFs", metavar="PATH")
parser.add_argument("-p", "--orfs", nargs=1, required=True, help="Protein identification result file", metavar="PATH")
parser.add_argument("-o", "--output", nargs=1, required=True, help="Filtered BLAST result outfile", metavar="PATH")

args = parser.parse_args()
print(args)

protRevStr="XXX_"
protContStr="CONT_"
pepTh=1
prtCSV = readFile(args.orfs[0], ',')
#prtCSV = readFile(prtFile, ',')
filteredProts=pepThresholding(prtCSV, pepTh, protRevStr, protContStr)
queryList=filteredProts['description']
filterBlast(args.blast[0], queryList, args.output[0])
#prtFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/pasa_transdecoder+fdr+th+grouping+prt.csv"
#blastFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/pasa_transdecoder_nonstar.csv"
#outBlast="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/pasa_transdecoder_nonstar_identified.csv"
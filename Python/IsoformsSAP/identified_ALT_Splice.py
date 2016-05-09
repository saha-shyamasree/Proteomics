##This code filters Identified ALT_SPLICE RNAs

import pandas as pd
import argparse

def readFile(filename, sep, header):
    if header==0:
        fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''], header=None)
    else:
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

def identifiedAlt_Splice(protList, altList):
    return protList[protList.isin(altList)]
    
protRevStr="XXX_"
protContStr="CONT_"
pepTh=1

parser = argparse.ArgumentParser(description='This script read output of UniProteinLocation.py and identify variations')
parser.add_argument("-p", "--protein", nargs=1, required=True, help="full path to the protein identification file", metavar="PATH")
parser.add_argument("-a", "--alt_splice", nargs=1, required=True, help="full path to the alt_splice_label_combinations.dat file", metavar="PATH")
args = parser.parse_args()

prots=readFile(args.protein[0],',',1)
alt_splice=readFile(args.alt_splice[0],'\t',0)

fProts=pepThresholding(prots, pepTh, protRevStr, protContStr)
protList=fProts['protein accession'].str.split('|').str.get(0)
alt_ids=alt_splice[1]

##This code to calculated number of canonical, isoform proteins in standard search. When we have multiple samples and
##MS data is splitted into multiple files according to samples, and MSGF search against standard database was cariied out
##sample specific way, this code will calculte total number of protein, canonical, isoform, TrEMBL and so on.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse
import glob
import pandas as pd
import numpy as np
import re

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;


parser = argparse.ArgumentParser(description='This python code counts unique number of proteins and peptides in a dataset')

parser.add_argument("-p", "--psm", nargs=1, required=True, help="full path of PSM and protein files folder", metavar="PATH")


args = parser.parse_args()
#print(args)
#onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
onlyfiles=glob.glob(args.psm[0]+"/*.standard+fdr+th+grouping+prt_filtered.csv")
#allProtein=readFile(args.fasta[0]+"TGEs.tsv", "\t")
#vcfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PITDB/Variations-proVCF/"
#proteinMat=pd.DataFrame(columns=('PAG ID','PAG score','protein accession','Pass Threshold (Protein)','description','group membership','distinct peptide sequences','ProteoGrouper:PDH score','unique peptides','razor peptides','cluster identifier'))
i=0
for f in onlyfiles:
    print("F:"+f)
    fBase=re.sub(re.escape(".standard+fdr+th+grouping+prt_filtered.csv"),"",f)
    #print("fBase"+fBase)
    sample=os.path.basename(fBase).split(".",maxsplit=1)[0]
    print("sample:"+str(sample))
    protDF=readFile(f, ',')
    pepDF=readFile(args.psm[0]+sample+".standard+fdr+th+grouping_filtered.csv", ',')
    if i==0:
        prot=protDF
        pep=pepDF
        i=i+1
    else:
        prot=prot.append(protDF)
        pep=pep.append(pepDF)

    
uniqueProtMat=prot.drop_duplicates(['protein accession'], keep='first')
uniqueProtMat.to_csv(args.psm[0]+"/uniqueProteinsStandard.csv")
uniquePepMat=pep.drop_duplicates(['Sequence'], keep='first')
uniquePepMat.to_csv(args.psm[0]+"/uniquePeptidesStandard.csv")
print("Total Protein:")
print(uniqueProtMat.shape)
print("Total Peptide:")
print(uniquePepMat.shape)
print("Reviewed Protein:")
print(uniqueProtMat.loc[(uniqueProtMat['protein accession'].str.contains("^sp"))].shape)
print("TrEMBL Protein:")
print(uniqueProtMat.loc[(uniqueProtMat['protein accession'].str.contains("^tr"))].shape)
print("Canonical Protein:")
print(uniqueProtMat.loc[(uniqueProtMat['protein accession'].str.contains("^sp")) & (~uniqueProtMat['protein accession'].str.contains("-"))].shape)
print("Known Isoform Protein:")
print(uniqueProtMat.loc[(uniqueProtMat['protein accession'].str.contains("^sp")) & (uniqueProtMat['protein accession'].str.contains("-"))].shape)

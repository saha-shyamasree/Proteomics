##This script read vcf file produced by peptideEvidenceIsoform.py and compute number of variations
##with peptide evidence, unique peptide and number of TGE observation with peptide evidence.

import pandas as pd
import argparse
import re
import os

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def tge_class_stats(vcf):
    ##this function counts peptides in respect to tge class
    print(vcf)
    
def stats(vcf, sample):
    #This function counts and print the results.
    doubleIds=vcf[vcf['ID'].astype(str).str.contains("\.2")]['QueryID']
    doubleEnded=vcf[vcf['QueryID'].isin(doubleIds)]
    doubleEndOne=doubleEnded[doubleEnded['ID'].astype(str).str.contains("\.1")]
    doubleEndTwo=doubleEnded[doubleEnded['ID'].astype(str).str.contains("\.2")]
    singleEnded=vcf[ ~vcf['QueryID'].isin(doubleIds)]
    singleEndedShortened=singleEnded[singleEnded['Type'].str.contains("shortened")]
    doubleEndedOneShortened=doubleEndOne[doubleEndOne['Type'].str.contains("shortened")]
    doubleEndedTwoShortened=doubleEndTwo[doubleEndTwo['Type'].str.contains("shortened")]
    doubleEndedShortenedIds=set(doubleEndedOneShortened['QueryID'].tolist()).intersection(doubleEndedTwoShortened['QueryID'].tolist())
    pepEvdTGE=len(vcf[vcf['PeptideCount'].astype(int)>0]['QueryID'].unique())
    pepEvdVar=vcf[vcf['PeptideCount'].astype(int)>0]['QueryID'].shape[0]
    unqPepEvdTGE=len(vcf[vcf['UniquePeptideCount'].astype(int)>0]['QueryID'].unique())
    unqPepEvdVar=vcf[vcf['UniquePeptideCount'].astype(int)>0]['QueryID'].shape[0]
    junctPep=vcf[(vcf.Evidence =="junction") | (vcf.Evidence =="juntion")]['QueryID'].unique().size
    #print("Sample\tTotal Isoforms\tSingle-sided\tDouble-sided\tShortened\tPeptideEvd(TGEObs)\tPeptideEvd(Variation)\tUniqueEvd(TGEObs)\tUniqueEvd(Variation)\n")
    print(sample+"\t"+str(len(vcf['QueryID'].unique()))+"\t"+str(singleEnded.shape[0])+"\t"+str(len(doubleEnded['QueryID'].unique()))+"\t"+str(len(doubleEndedShortenedIds)+singleEndedShortened.shape[0])+"\t"+str(pepEvdTGE)+"\t"+str(pepEvdVar)+"\t"+str(unqPepEvdTGE)+"\t"+str(unqPepEvdVar)+"\t"+str(junctPep))

def main(vcfFile): #, annotFile
    ##This is main function
    #annot=readFile(annotFile,',')
    #isoforms=annotations.loc[annotations['Class'].str.contains("prime_")]
    #isoIds=isoforms['ORF Id'].str.extract("([^ ]*)").tolist()
    sample=os.path.basename(vcfFile).replace(".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf","")
    #print("sample:"+sample)
    vcf=readFile(vcfFile,'\t')
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Evidence','Score'])
    #SubjectId=sp|Q92466-5|DDB2_HUMAN;QueryId=asmbl_10851|m.134573;QueryLength=244;QueryStart=1;QueryEnd=237;SubjectLength=244;SubjectStart=1;SubjectEnd=237;Type=3prime_alternative_C-terminus_partial;QPOS=237
    
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    #print(vcf.QueryID[0:5])
    vcf.QueryLength=vcf.QueryLength.str.replace('QueryLength=','')
    vcf.QueryStart=vcf.QueryStart.str.replace('QueryStart=','')
    vcf.QueryEnd=vcf.QueryEnd.str.replace('QueryEnd=','')
    vcf.SubjectLength=vcf.SubjectLength.str.replace('SubjectLength=','')
    vcf.SubjectStart=vcf.SubjectStart.str.replace('SubjectStart=','')
    vcf.SubjectEnd=vcf.SubjectEnd.str.replace('SubjectEnd=','')
    vcf.Type=vcf.Type.str.replace('Type=','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS=','')
    vcf.PeptideCount=vcf.PeptideCount.str.replace('PeptideCount=','')
    vcf.UniquePeptideCount=vcf.UniquePeptideCount.str.replace('UniquePeptideCount=','')
    vcf.Peptides=vcf.Peptides.str.replace('Peptides=','')
    vcf.Evidence=vcf.Evidence.str.replace('Evidence=','')
    vcf.Score=vcf.Score.str.replace('Score=','')
    stats(vcf, sample)
parser = argparse.ArgumentParser(description='This script read output of peptideEvidenceIsoforms.py and count variations')
parser.add_argument("-v", "--vcf", nargs=1, required=True, help="full path to the isoform peptide vcf file", metavar="PATH")
#parser.add_argument("-a", "--annot", nargs=1, required=True, help="full path to the annotation file", metavar="PATH")
args = parser.parse_args()
main(args.vcf[0]) #, args.annot[0]

#vcfFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Annotation/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf"
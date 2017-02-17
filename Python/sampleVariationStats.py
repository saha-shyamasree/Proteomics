## This code reads Olivers ORFs and transcript fasta files, and protein and peptide identification files.
## Filters ORFs and transcript fasta files using fastaFileFiltering.py. It also filters the peptide csv.

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import subprocess
import re
import argparse
import pandas as pd

## Read Identified Proteins

def readIdentifiedProteinPeptide(filename,sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj

def readFastaFile(fastaFile):
    ##This function reads a fasta file using SeqIO, convert entries into record object and put in dataframe.
    fastaHandle = open(fastaFile, "rU")
    records = list(SeqIO.parse(fastaHandle, "fasta"))
    #print("Records:"+str(len(records)))
    fastaDF=pd.DataFrame(columns=("ORF Id","Sequence"))
    for r in records:
        fastaDF=fastaDF.append({'ORF Id':r.description,'Sequence':str(r.seq)}, ignore_index=True)
        #print("seqList:"+str(len(seqList)))
    return fastaDF

def readFastaSeq(fastaFile):
    fastaHandle = open(fastaFile, "rU")
    records = list(SeqIO.parse(fastaHandle, "fasta"))
    seqList=[str(r.seq) for r in records]
    return list(set(seqList))

def vcfVarReader(filename):
    vcf=readIdentifiedProteinPeptide(filename, '\t')
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Score'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
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
    vcf.Score=vcf.Score.str.replace('Score=','')
    return vcf


def separateAltSplicingFromINDELs(vcfObj):
    vcfObj['RefCount']=vcfObj['REF'].str.len()
    vcfObj['ALTCount']=vcfObj['ALT'].str.len()
    isoVCF=vcfObj[(vcfObj['RefCount']>9) | (vcfObj['ALTCount']>9)]
    uniqIsoIds=isoVCF['QueryID'].unique().tolist()
    vcfObj.loc[(vcfObj['RefCount']>9) | (vcfObj['ALTCount']>9),'Type']="ALT_SPLICE"
    return uniqIsoIds

def printVarCounts(sample,unqVCFMat):
    ##this function prints SAPs/ALT/INDELs and peptide evidence count for those.
    #print(unqVCFMat.shape)
    ##SAP Counts
    sap=str(unqVCFMat[unqVCFMat['Type']=="SAP"].shape[0])
    psap=str(unqVCFMat[(unqVCFMat['Type']=="SAP") & (unqVCFMat['FILTER']=="PASS")].shape[0])
    ssap=str(unqVCFMat[unqVCFMat['Type']=="SSAP"].shape[0])
    pssap=str(unqVCFMat[(unqVCFMat['Type']=="SSAP") & (unqVCFMat['FILTER']=="PASS")].shape[0])
    sorf=str(len(unqVCFMat.loc[unqVCFMat['Type']=="SAP"]['QueryID'].unique()))
    ssorf=str(len(unqVCFMat.loc[unqVCFMat['Type']=="SSAP"]['QueryID'].unique()))
    
    ##ALT counts
    alt=str(unqVCFMat[unqVCFMat['Type']=="ALT"].shape[0])
    palt=str(unqVCFMat[(unqVCFMat['Type']=="ALT") & (unqVCFMat['FILTER']=="PASS")].shape[0])
    salt=str(unqVCFMat[unqVCFMat['Type']=="SALT"].shape[0])
    psalt=str(unqVCFMat[(unqVCFMat['Type']=="SALT") & (unqVCFMat['FILTER']=="PASS")].shape[0])
    aorf=str(len(unqVCFMat.loc[unqVCFMat['Type']=="ALT"]['QueryID'].unique()))
    saorf=str(len(unqVCFMat.loc[unqVCFMat['Type']=="SALT"]['QueryID'].unique()))
    #print(unqVCFMat[unqVCFMat['Type']=="ALT"].shape)
    ##INS counts
    ins=str(unqVCFMat[unqVCFMat['Type']=="INS"].shape[0])
    pins=str(unqVCFMat[(unqVCFMat['Type']=="INS") & (unqVCFMat['FILTER']=="PASS")].shape[0])
    iorf=str(len(unqVCFMat.loc[unqVCFMat['Type']=="INS"]['QueryID'].unique()))
    
    ##DEL counts
    dele=str(unqVCFMat[unqVCFMat['Type']=="DEL"].shape[0])
    pdele=str(unqVCFMat[(unqVCFMat['Type']=="DEL") & (unqVCFMat['FILTER']=="PASS")].shape[0])
    dorf=str(len(unqVCFMat.loc[unqVCFMat['Type']=="DEL"]['QueryID'].unique()))
    
    printVar=sample+"\t"+ssorf+"\t"+ssap+"\t"+pssap+"\t"+sorf+"\t"+sap+"\t"+psap+"\t"+saorf+"\t"+salt+"\t"+psalt+"\t"+aorf+"\t"+alt+"\t"+palt+"\t"+iorf+"\t"+ins+"\t"+pins+"\t"+dorf+"\t"+dele+"\t"+pdele
    print(printVar)

def mainVar(args):
    sample=re.sub(re.escape(".assemblies.fasta.transdecoder.pep_pepEvd.vcf"),"",os.path.basename(args.vcf))
    ##Read annotation file to see how many are known, novel isoform, variations
    protAnnot=readIdentifiedProteinPeptide(args.annot,',')
    ##reads mapping zone variation vcf File to determine which of the ORFs should be considered as isoform and not known variation
    vcfVar=vcfVarReader(args.vcf)
    
    vcfIsoIds=separateAltSplicingFromINDELs(vcfVar)
    protAnnot['mRNA']=protAnnot['ORF Id'].str.extract("([^ ]+)",expand=False)
    
    protAnnot['NewClass']=protAnnot['Class']
    protAnnot.loc[protAnnot['mRNA'].isin(vcfIsoIds),'NewClass']="ALT_SPLICE"
    #print("Total vcf ids")
    #print(len(vcfVar['QueryID'].unique()))
    #print("vcf id in annot")
    #print(len(vcfVar.loc[vcfVar['QueryID'].isin(protAnnot['mRNA'])]['QueryID'].unique()))
    ##Reading identified ORFs and storing that in panda dataframe object
    protDF=readFastaFile(args.ORFsIdent)
    protMerged=pd.merge(protAnnot,protDF, on='ORF Id', how='outer')
    #protMerged=pd.merge(protMerged,vcfVarCount, on='mRNA', how='outer')
    uniqueProtMat=protMerged.drop_duplicates(['Protein ID', 'Class', 'Variation', 'Species','Protein Name', 'Gene Name', 'Protein description', 'Source', "NewClass",'Sequence'], keep='first')
    unqIds=uniqueProtMat['mRNA']
    uniqVCFVarIso=vcfVar[vcfVar['QueryID'].isin(unqIds)]
    #uniqVCFVarIso.to_csv("uniqVCFVarIso.csv",sep="\t")
    ##remove variations that are longer than 9AAs
    #uniqVCFVar=uniqVCFVarIso[uniqVCFVarIso['Type']!='ALT_SPLICE']
    #uniqVCFVar.to_csv("uniqVCFVar.csv",sep="\t")
    printVarCounts(sample,uniqVCFVarIso)
    

parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("--ORFsIdent", required=True, help="Identified ORFs fasta file name")
parser.add_argument("--annot", required=True, help="Protein csv out file name")
parser.add_argument("--vcf", required=True, help="Peptide csv out file name")

args = parser.parse_args()
mainVar(args)

    
## Test Command
##python sampleVariationStats.py --ORFsIdent /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta --annot /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Annotation/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Variations-proVCF/human_adeno.assemblies.fasta.transdecoder.pep_pepEvd.vcf"
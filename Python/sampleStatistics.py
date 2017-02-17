## This code reads ORFs and transcript fasta files, and protein and peptide identification files.
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
    
def identifiedProteinList(prots, columnName):
    identifiedProts=prots[columnName].tolist()
    return identifiedProts
    
def extractProtIdsFromDescription(prots, sep):
    ##prots is a list of strings, generally the description column of protein identification file
    protIds=[t.partition(sep)[0] for t in prots]
    return protIds

def writeResult(filteredObj, filename):
    filteredObj.to_csv(filename, index=False)
def vcfReader(filename):
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
    return uniqIsoIds

def main(args):
    sample=re.sub(re.escape("+fdr+th+grouping_filtered.csv"),"",os.path.basename(args.peptides))
    #print("sample:"+sample)
    proc1=subprocess.check_output(["grep \"TITLE\" "+args.ms+" | wc -l"], shell=True)
    totalSpec=proc1.decode("utf-8").strip()
    protColName="description"
    pepColName="proteinacc_start_stop_pre_post_;"
    pepObj=readIdentifiedProteinPeptide(args.peptides,",")
    
    identSpec=len(pepObj['Spectrum Title'].unique())
    proc2 = subprocess.check_output(["grep '>' "+args.transcripts+" | wc -l"], shell=True)
    allTranscript=proc2.decode("utf-8").strip()
    #print("allTrans:"+allTranscript)
    proc3 = subprocess.check_output(["grep '>' "+args.transcriptsIdent+" | wc -l"], shell=True)
    identifiedTrans=proc3.decode('utf-8').strip()
    all_ORFs=readFastaSeq(args.ORFs)#readFastaFile(args.ORFs)
    proc4=subprocess.check_output(["awk -F , '{if (NR!=1) {print $1}}' "+args.proteins+"| sort | uniq | wc -l"], shell=True)
    pagCount=proc4.decode('utf-8').strip()
    #ident_ORFs=readFastaSeq(args.ORFsIdent)
    ##Read annotation file to see how many are known, novel isoform, variations
    protAnnot=readIdentifiedProteinPeptide(args.annot,',')
    ##reads vcf File to determine which of the ORFs should be considered as isoform and not known variation
    vcf=vcfReader(args.vcf)
    vcfIsoIds=separateAltSplicingFromINDELs(vcf)
    protAnnot['mRNA']=protAnnot['ORF Id'].str.extract("([^ ]+)",expand=False)
    protAnnot['NewClass']=protAnnot['Class']
    protAnnot.loc[protAnnot['mRNA'].isin(vcfIsoIds),'NewClass']="ALT_SPLICE"
    protDF=readFastaFile(args.ORFsIdent)
    protMerged=pd.merge(protAnnot,protDF, on='ORF Id', how='outer')
    uniqueProtMat=protMerged.drop_duplicates(['Protein ID', 'Class', 'Variation', 'Species','Protein Name', 'Gene Name', 'Protein description', 'Source', "NewClass",'Sequence'], keep='first')
    identORFs=uniqueProtMat.shape[0]
    canonical=uniqueProtMat.loc[(uniqueProtMat['Class']=='known') & (uniqueProtMat['Source']=='sp') & ~(uniqueProtMat['Protein ID'].str.contains('-'))].shape[0]
    trembl=uniqueProtMat.loc[(uniqueProtMat['Class']=='known') & (uniqueProtMat['Source']=='tr')].shape[0]
    knownIso=uniqueProtMat.loc[(uniqueProtMat['Class']=='known') & (uniqueProtMat['Source']=='sp') & (uniqueProtMat['Protein ID'].str.contains('-'))].shape[0]
    variation=uniqueProtMat.loc[uniqueProtMat['NewClass']=="known variation"].shape[0]
    novel=uniqueProtMat.loc[uniqueProtMat['NewClass']=="novel"].shape[0]
    novelIso=uniqueProtMat.loc[(uniqueProtMat['NewClass']!='known') & (uniqueProtMat['NewClass']!='novel') & (uniqueProtMat['NewClass']!="known variation")].shape[0]
    
    ##Print Info
    printStr=sample+"\t"+str(totalSpec)+"\t"+str(identSpec)+"\t"+"READ_COUNT"+"\t"+str(allTranscript)+"\t"+str(identifiedTrans)+"\t"+str(len(all_ORFs))+"\t"+str(len(pepObj['Sequence'].unique().tolist()))+"\t"+pagCount+"\t"+str(identORFs)+"\t"+str(canonical)+"\t"+str(knownIso)+"\t"+str(trembl)+"\t"+str(variation)+"\t"+str(novelIso)+"\t"+str(novel)
    print(printStr)
    

parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("--ORFs", required=True, help="ORFs fasta file name")
parser.add_argument("--transcripts", required=True, help="Transcript fasta file name")
parser.add_argument("--proteins", required=True, help="mzIdentML-lib protein export file")
parser.add_argument("--peptides", required=True, help="mzIdentML-lib PSM export file")
parser.add_argument("--ORFsIdent", required=True, help="Identified ORFs fasta file name")
parser.add_argument("--transcriptsIdent", required=True, help="Identified Transcripts fasta file name")
parser.add_argument("--annot", required=True, help="Protein csv out file name")
parser.add_argument("--vcf", required=True, help="Peptide csv out file name")
parser.add_argument("--ms", required=True, help="Peptide csv out file name")

args = parser.parse_args()
main(args)

    
## Test Command
## python PIT-DBData_processing.py --ORFs D:\data\Oliver\ORFs\G10.assemblies.fasta.transdecoder.pep --transcripts D:\data\Oliver\transcripts\G10.assemblies.fasta --proteins D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.pep+fdr+th+grouping+prt.csvDBIds.csv --peptides D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.pep+fdr+th+grouping.csv
##--ORFsOutFile D:\data\Oliver\ORFs\G10.assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile D:\data\Oliver\transcripts\G10.assemblies.fasta.identified.fasta --peptideOutFile D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.pep+fdr+th+grouping_filtered.csv
##python sampleStatistics.py --transcripts Statistics.py --transcripts "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.assemblies.fasta" --peptides "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping_filtered.csv" --transcriptsIden "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.assemblies.fasta.identified.fasta" --ms "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/MS/DM_from_raw.mgf"
##python sampleStatistics.py --ORFs /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep --transcripts "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.assemblies.fasta" --peptides "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping_filtered.csv" --transcriptsIden "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.assemblies.fasta.identified.fasta" --ms "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/MS/DM_from_raw.mgf" --proteins "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/PSMs-Peptides-ORFs/human_adeno+fdr+th+grouping+prt_filtered.csv" --ORFsIdent /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta
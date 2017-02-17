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
from itertools import chain

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

def vcfIsoReader(filename):
    vcf=readIdentifiedProteinPeptide(filename, '\t')
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Evidence','Score'])
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
    vcf.Evidence=vcf.Evidence.str.replace('Evidence=','')
    vcf.Score=vcf.Score.str.replace('Score=','')
    return vcf

def separateAltSplicingFromINDELs(vcfObj):
    isoVCF=vcfObj[(vcfObj['RefCount']>9) | (vcfObj['ALTCount']>9)]
    uniqIsoIds=isoVCF['QueryID'].unique().tolist()
    vcfObj.loc[(vcfObj['RefCount']>9) | (vcfObj['ALTCount']>9),'TYPE']="ALT_SPLICE"
    return uniqIsoIds

def varToIso(vcfObj):
    isoVCF=vcfObj[(vcfObj['RefCount']>9) | (vcfObj['ALTCount']>9)]
    return isoVCF


def printIsoCounts(sample, vcfIso, novelIso, isoVar, isoClasses):
    ##this function prints number of TGEs in each isoform classes
    countDF=pd.DataFrame(columns=['Class','Count','Peptide'])
    #countDF['Class']=isoClasses
    #unqNovelIso=novelIso['mRNA']
    #alt_spliced=novelIso.loc[novelIso['NewClass']=='ALT_SPLICE']['mRNA']
    #print(alt_spliced)
    #altVCFIso=vcfIso.loc[vcfIso['mRNA'].isin(alt_spliced)]
    #print(altVCFIso)
    #print("Total vcf Iso:"+str(len(vcfIso['mRNA'].unique())))
    for i in range(0,len(isoClasses)):
        iso=isoClasses[i]
        count=int(novelIso[novelIso['NewClass']==iso].shape[0])
        novelIsoIds=novelIso.loc[novelIso['NewClass']==iso]['mRNA']
        #print(novelIsoIds)
        if iso=="ALT_SPLICE":
            #in initial classification these lonf ALT/INDELs events were considered as variation, but because they are
            #longer than 9AAs, we change thir class from known variation to alt_splice [note some isoform may also have
            #these events but they will still be classified according to their boundary alterations]. Hence peptide evidence
            #for these TGEs can be found in variation vcf and not in the isoform vcf
            peptide=int(len(isoVar.loc[(isoVar['QueryID'].isin(novelIsoIds)) & (isoVar['FILTER']=='PASS')]['QueryID'].unique()))
            #print(peptide)
            #print(isoVar.loc[(isoVar['QueryID'].isin(novelIsoIds)) & (isoVar['FILTER']=='PASS')]['QueryID'].unique())
        else:
            peptide=int(len(vcfIso.loc[(vcfIso['mRNA'].isin(novelIsoIds)) & (vcfIso['FILTER']=='PASS')]['mRNA'].unique()))
        
        countDF.loc[i]=[iso,count,peptide]
        count=0
        peptide=0
    #print(countDF.to_csv())
    #print("\t".join(countDF['Class'].tolist()).replace("prime_",""))
    ##ADD sample name here
    countDF.Count=countDF.Count.astype(int)
    countDF.Peptide=countDF.Peptide.astype(int)
    #print(countDF[['Count','Peptide']])
    val=countDF[['Count','Peptide']].values.tolist()
    unlistVal=str(list(chain(*val)))
    #print(str(unlistVal))
    #uVal=[str(v) for v in unlistVal]
    printVar=sample+"\t"+unlistVal
    printVar=printVar.replace(",","\t")
    printVar=printVar.replace("[","")
    printVar=printVar.replace("]","")
    print(printVar)

def mainIsoform(args):
    sample=re.sub(re.escape(".assemblies.fasta.transdecoder.pep_pepEvd.vcf"),"",os.path.basename(args.vcf))
    ##Read annotation file to see how many are known, novel isoform, variations
    protAnnot=readIdentifiedProteinPeptide(args.annot,',')
    ##reads mapping zone variation vcf File to determine which of the ORFs should be considered as isoform and not known variation
    vcfVar=vcfVarReader(args.vcf)
    vcfVar['RefCount']=vcfVar['REF'].str.len()
    vcfVar['ALTCount']=vcfVar['ALT'].str.len()
    vcfIsoIds=separateAltSplicingFromINDELs(vcfVar)
    #print(vcfIsoIds)
    protAnnot['mRNA']=protAnnot['ORF Id'].str.extract("([^ ]+)",expand=False)
    protAnnot['NewClass']=protAnnot['Class']
    #protAnnot['Evidence']='No'
    protAnnot.loc[(protAnnot['mRNA'].isin(vcfIsoIds)) & (protAnnot['Class']=='known variation'),'NewClass']="ALT_SPLICE"
    isoVar=varToIso(vcfVar)
    ##Read isoform vcf file to check how many of these has peptide evidence
    
    vcfIso=vcfIsoReader(args.vcfIso)
    vcfIso['mRNA']=vcfIso['QueryID'].str.extract("([^ ]+)")    
    ##Reading identified ORFs and storing that in panda dataframe object
    protDF=readFastaFile(args.ORFsIdent)
    protMerged=pd.merge(protAnnot,protDF, on='ORF Id', how='outer')
    #protMerged=pd.merge(protMerged,vcfVarCount, on='mRNA', how='outer')
    uniqueProtMat=protMerged.drop_duplicates(['Protein ID', 'Class', 'Variation', 'Species','Protein Name', 'Gene Name', 'Protein description', 'Source', "NewClass",'Sequence'], keep='first')
    unqIds=uniqueProtMat['mRNA']
    isoClasses=['5prime_shortened','3prime_shortened','5prime_extended','3prime_extended','5prime_alternative','3prime_alternative','5prime_shortened_3prime_shortened','5prime_shortened_3prime_extended','5prime_shortened_3prime_alternative','5prime_extended_3prime_shortened','5prime_extended_3prime_extended','5prime_extended_3prime_alternative','5prime_alternative_3prime_shortened','5prime_alternative_3prime_extended','5prime_alternative_3prime_alternative','ALT_SPLICE']
    novelIso=uniqueProtMat[(uniqueProtMat['NewClass']!='known') & (uniqueProtMat['NewClass']!='novel') & (uniqueProtMat['NewClass']!="known variation")]
    printIsoCounts(sample, vcfIso, novelIso, isoVar, isoClasses)

    
    

parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("--ORFsIdent", required=True, help="Identified ORFs fasta file name")
parser.add_argument("--annot", required=True, help="Protein csv out file name")
parser.add_argument("--vcf", required=True, help="Peptide csv out file name")
parser.add_argument("--vcfIso", required=True, help="Peptide csv out file name")

args = parser.parse_args()

mainIsoform(args)
    
## Test Command
##python sampleIsoformStats.py --ORFsIdent /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.assemblies.fasta.transdecoder.pep.identified.fasta --annot /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Annotation/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Variations-proVCF/human_adeno.assemblies.fasta.transdecoder.pep_pepEvd.vcf" --vcfIso "/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Annotation/human_adeno_mydb_pasa.assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf"
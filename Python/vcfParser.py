###vcf parser
import pandas as pd
import csv
import re
import argparse

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def vcfParser(vcfFileName):
    ##read vcf file
    vcf=readFile(vcfFileName,'\t')

    ## split the INFO column and add them to main data frame
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Score'])
    #print(vcf.QueryID[0:5])
    
    info.Alignment=info.Alignment.str.replace('Alignment=\[','')
    info.Alignment=info.Alignment.str.replace('\]','')
    ##Split Alignemnt column
    alignment=pd.DataFrame(info.Alignment.str.split(':').tolist(),columns=['QueryLength','QueryStart','QueryEnd','SubjectLength','SubjectStart','SubjectEnd'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info['SubjectID'])
    vcf=vcf.join(info['QueryID'])
    vcf=vcf.join(alignment)
    vcf=vcf.join(info['Type'])
    vcf=vcf.join(info['QPOS'])
    vcf=vcf.join(info['PeptideCount'])
    vcf=vcf.join(info['UniquePeptideCount'])
    vcf=vcf.join(info['Peptides'])
    vcf=vcf.join(info['Score'])
    
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    vcf.QueryLength=vcf.QueryLength.str.replace('QueryLength=','')
    vcf.QueryStart=vcf.QueryStart.str.replace('QueryStart=','')
    vcf.QueryEnd=vcf.QueryEnd.str.replace('QueryEnd=','')
    vcf.SubjectLength=vcf.SubjectLength.str.replace('SubjectLength=','')
    vcf.SubjectStart=vcf.SubjectStart.str.replace('SubjectStart=','')
    vcf.SubjectEnd=vcf.SubjectEnd.str.replace('SubjectEnd=','')    
    vcf.Type=vcf.Type.str.replace('Type:','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS:','')
    vcf.PeptideCount=vcf.PeptideCount.str.replace('PeptideCount:','')
    vcf.UniquePeptideCount=vcf.UniquePeptideCount.str.replace('UniquePeptideCount:','')
    vcf.Peptides=vcf.Peptides.str.replace('Peptides:','')
    vcf.Score=vcf.Score.str.replace('Score:','')
    return vcf

parser = argparse.ArgumentParser(description='This script read output of UniProteinLocation.py and identify variations')
parser.add_argument("-v", "--vcf", nargs=1, required=True, help="full path to the vcf file", metavar="PATH")

args = parser.parse_args()
print(args.vcf[0])

vcfObj=vcfParser(args.vcf[0])
print(vcfObj.to_csv())
print(vcfObj.columns.values)
#return vcfObj

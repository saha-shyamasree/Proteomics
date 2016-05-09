## This code counts variations per type, i.e. SAP/SAAV, ALT, INDELs. Ir reads 2 vcf files, one for identified ORFs
## another for variations from identified ORFs with peptide evidence.

import pandas as pd
import argparse

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def stats(vcf,flag):
    if flag==1:
        info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS'])
    elif flag==2:
        info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Score'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    vcf.Alignment=vcf.Alignment.str.replace('Alignment=\[','')
    vcf.Alignment=vcf.Alignment.str.replace('\]','')
    vcf.Type=vcf.Type.str.replace('Type:','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS:','')
    vcfGrouped=vcf.groupby('Type')
    for name, variations in vcfGrouped:
        print(name+":"+str(len(variations)))

    
parser = argparse.ArgumentParser(description='This code counts variations per type, i.e. SAP/SAAV, ALT, INDELs. Ir reads 2 vcf files, one for identified ORFs and another for variations from identified ORFs with peptide evidence.')

parser.add_argument("-v", "--vcf", nargs=1, required=True, help="full path to the all vcf file", metavar="PATH")
parser.add_argument("-o", "--vcfident", nargs=1, required=True, help="full path to the peptide evidence vcf file", metavar="PATH")


args = parser.parse_args()

allVcf=readFile(args.vcf[0],'\t')
stats(allVcf,1)

identVcf=readFile(args.vcfident[0],'\t')
stats(identVcf,2)
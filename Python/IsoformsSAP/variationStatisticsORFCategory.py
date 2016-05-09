## This code counts variations per type, i.e. SAP/SAAV, ALT, INDELs for protein categories i.e. known protein variation and isoform variation.
##It reads 2 vcf files, one for identified ORFs another for variations from identified ORFs with peptide evidence.

import pandas as pd
import argparse


class AminoAcidVariation:
    """This clas holds aminoacid variations from blast alignment"""
    def __init__(self,subjectId, queryId, start, end, Type, vid):
        self.subjectId=subjectId
        self.queryId=queryId
        self.pos=start
        self.qpos=qpos
        self.ref=refSeq
        self.alt=altSeq
        self.type=Type
        self.chr=chrm
        self.vid=vid
        self.alingmentInfo=info
    def getSubjectId(self):
        return self.subjectId
    def getQueryId(self):
        return self.queryId
    def getPos(self):
        return self.pos

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def statsORFCategoryWise(vcf,flag, protList):
    print("in function1")
    ##this function identifies number of variations cominng from the known protein variation and the isoform with variation
    if flag==1:
        info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS'])
    elif flag==2:
        info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Score'])
        protList=protList.str.replace('\s+',' ')
        ### following line is necessary to accomodate the fact that the ORFs ids produced by transdecoder has spaces and the MSGF identification consider the [^\s]+ as the id of the ORF. To match theses ids we remove everything after the first space.
        queryIDs=pd.DataFrame(protList.str.split(' ').tolist(),columns=['QueryID','GeneID','ORF','GeneID2','QueryID2','Type','Length','Strand','Location'])
        protList=queryIDs['QueryID']
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    vcf.Alignment=vcf.Alignment.str.replace('Alignment=\[','')
    vcf.Alignment=vcf.Alignment.str.replace('\]','')
    vcf.Type=vcf.Type.str.replace('Type:','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS:','')
    
    ##If the variation is coming from one of the proteins from the protList, we count it.
    identIdx=vcf.QueryID.isin(protList)
    idx=identIdx[identIdx==True].index.tolist()
    vcfIdentified= vcf.iloc[idx,:]
    
    vcfGrouped=vcfIdentified.groupby('Type')
    for name, variations in vcfGrouped:
        print(name+":"+str(len(variations)))
    
    print("in function2")

  
parser = argparse.ArgumentParser(description='This code counts variations per type, i.e. SAP/SAAV, ALT, INDELs. Ir reads 2 vcf files, one for identified ORFs and another for variations from identified ORFs with peptide evidence.')

parser.add_argument("-v", "--vcf", nargs=1, required=True, help="full path to the all vcf file", metavar="PATH")
parser.add_argument("-o", "--vcfident", nargs=1, required=True, help="full path to the peptide evidence vcf file", metavar="PATH")
parser.add_argument("-p", "--protein", nargs=1, required=True, help="full path to the protein file", metavar="PATH")

args = parser.parse_args()

proteinList=readFile(args.protein[0],',')

allVcf=readFile(args.vcf[0],'\t')
statsORFCategoryWise(allVcf,1,proteinList['ORF Id'])


identVcf=readFile(args.vcfident[0],'\t')

statsORFCategoryWise(identVcf,2,proteinList['ORF Id'])
print("peptide evidence")
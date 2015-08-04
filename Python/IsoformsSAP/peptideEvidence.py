## This code reads the SAPs/Alt/INDELs vcf generated from IdentifyProteinIsoformSAP.py, list of known proteins,
## known proteins with SAPs, Isoforms, and Isoforms with SAPs, and the peptide identification result file from
## MSGF+ after post-processing using mzIdenML-lib. This scripts generates a vcf file containing SAPs/ALTs/INDELs
## with peptide evidence.
import pandas as pd
import csv
import re
from AminoAcidVariation import AminoAcidVariation

def checkVariationPeptideCoverage(protId, loc, sequence, altSeq):
    ## Check where this peptide is covers the location of the variation.
    ## If yes, compare the peptide sequence and the variation sequence to make sure the peptide does
    ## support the variation.
    ## Split the protId, get start and end and compare with the loc. Check the sequence, compare with altSeq.
    protList=protId.split(';')
    
def findPeptide(peptideGrouped, ORFId, loc, altSeq):
    PSMsList={}
    for seq, query in peptideGrouped:
        # here query holds all the PSMs for the peptide Sequence. 'proteinacc_start_stop_pre_post_;' column
        # lists the proteins it might have come from and its location with each of the proteins.
        # Accessing the first row of the group is enough as same sequence will have same protein list.
        protList=query.iloc[0]['proteinacc_start_stop_pre_post_;']
        if ORFId in protList:
            ## this means the ORF in question was identified and containes SAP/ALT/INDELs
            ## This prptide/Sequence is an evidence of the ORF in interest. Hence, put these PSMs together.
            ## Essentially this is grouping PSMs according to proteins. As protein file tells us which
            ## peptides belong to the sub-members of a PAG
            ## create a key value structure, where the orfsid is the key
            checkVariationPeptideCoverage(protList, loc, sequence, altSeq)
            PSMsList[ORFId]=query
            
  
def findPeptideEvidence(vcf, PSMs, newVcfFileName):
    ## vcf file fields: Chr, POS within Protein	ID, REF, ALT, INFO(SubjectId=P09417-2;QueryId=Dataset_A_asmbl_41426_ORF20_Frame_3_84-446;Alignment=[QueryLength=121:QueryStart=1:QueryEnd=116:SubjectLength=213:SubjectStart=1:SubjectEnd=116];Type:SSAP)
    ## split the INFO column and add them to main data frame
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryID=','')
    vcf.Alignment=vcf.Alignment.str.replace('Alignment=[','')
    vcf.Alignment=vcf.Alignment.str.replace(']','')
    ##groups vcf entries by ORF/Query ids.
    vcfGrouped=vcf.groupby('QueryID')
    
    ##in the same way group the PSMs according to the prptide sequence.
    peptideGrouped=PSMs.groupby('Sequence')
    # each of these group represents all the SAPs/ALTs/INDELs for a ORF/Query
    for name, query in vcfGrouped:
        ## check whether this ORF has been identified and whether this SAP/ALT/INDEL event has a peptide evidence.
        ## the ORF id contained ',' which was replaced by ';' in contigStat.pl as this code produce comma separated file. Later
        ## on IdentifyProeinIsoformSAP.py replaced ';' by '&' for same cause. ORFs with multiple parent
        ## transcripts have '&'. Whereas, the same ORFs in the PSMs list contain ','.
        pattern=re.compile('&')
        nameComma=pattern.sub(',',name)
        ## find peptide evidence of this ORF
        PSMsOfORF=findPeptide(peptideGrouped, ORFId)
        
        ## Match 
    
def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep)
    return fileDFObj;

def main(PSMFileName, vcfFileName, newVcfFile):
    ##read peptide identification file
    PSMs=readFile(PSMFileName, ',')
    ##read vcf file
    vcf=readFile(vcfFileName,'\t')
    ##for each entry in the vcf try find peptides that overlaps the vcf entry location. For deletion event, it might
    ##be little tricky. Not finding any peptide for deletion event is a good sign but does not confirms the deletion.
    findPeptideEvidence(vcf, PSMs, newVcfFileName)
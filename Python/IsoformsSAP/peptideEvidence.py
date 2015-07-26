## This code reads the SAPs/Alt/INDELs vcf generated from IdentifyProteinIsoformSAP.py, list of known proteins,
## known proteins with SAPs, Isoforms, and Isoforms with SAPs, and the peptide identification result file from
## MSGF+ after post-processing using mzIdenML-lib. This scripts generates a vcf file containing SAPs/ALTs/INDELs
## with peptide evidence.
import pandas as pd
import csv
from AminoAcidVariation import AminoAcidVariation

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
        ## check whether this protein has been identified and whether this SAP/ALT/INDEL event has a peptide evidence.
        ## find peptide evidence of this ORF
        PSMsOfORF=findPeptide(PSMs, ORFId)
        
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
## This code reads the SAPs/Alt/INDELs vcf generated from IdentifyProteinIsoformSAP.py, list of known proteins,
## known proteins with SAPs, Isoforms, and Isoforms with SAPs, and the peptide identification result file from
## MSGF+ after post-processing using mzIdenML-lib. This scripts generates a vcf file containing SAPs/ALTs/INDELs
## with peptide evidence.
##             [Many to Many]
##   Protein <#############> PSM 
##     #      #           #  ^  
##     #       #         #   #
##     #        #       #    #
##     #         #     #     #
##     #          #   #      #
##     #           # #       #
## [One#to Many]    #        # [Many to Many]
##     #           # #       #
##     #          #   #      #
##     #         #     #     #
##     #        #       #    #
##     #       #         #   #
##     V      #           #  V
##  Variation<#############> Petide
##           [Many to Many]


import pandas as pd
import csv
import re
from AminoAcidVariation import AminoAcidVariation

def checkVariationPeptideCoverage(ORFId, PSMsOfORF, variation):
    prevPeptide=''
    prevFound=0
    PSMEvidence=pd.DataFrame();
    for i in 0:len(PSMsOfORF):
        protAcc=PSMsOfORF.iloc[i]['proteinacc_start_stop_pre_post_;']
        if prevPeptide!=protAcc:
            ## Check whether this PSM covers the location of the variation.
            ## If yes, compare the peptide sequence and the variation sequence to make sure the peptide does
            ## support the variation.
            ## Split the protId, get start and end and compare with the loc. Check the sequence, compare with altSeq.
            prevPeptide=protAcc;
            ## prevFound flag is reset here. Because prevPeptide is not same as current one, therefore the PSM has not been
            ## added to the PSMRvidence yet.
            prevFound=0
            protList=protAcc.split(';')
            proacc_start_stop_pre_post=[s for s in protList if ORFId in s] #this should always have one element.
            if len(proacc_start_stop_pre_post)==1:
                ##this should always be the case.
                pattern=re.compile('(.*)?-(\d+)_(\d+)_(\d+)_([A-Z]|-)_([A-Z]|-)')
                match=pattern.finditer(proacc_start_stop_pre_post[0])
                count=0
                for m in match:
                    ##this loop should run only once.
                    proacc_start_stop_pre_postList=m
                    if count>0:
                        ##ERROR
                        print("ERROR: Match object should have single value.")
                        count=count+1;
                        break;
                    count=count+1;
                ## this is what should happen
                if count==1:
                    ## Check whether location of the SAP/ALT/INDEL is covered by this peptide.
                    ## For SAPs its straight forward. 
                    ## It might happen that for ALT/INDELs the peptide covers the event partially.
                    if variation['Type']=='SAP' or variation['Type']=='SSAP':
                        if int(proacc_start_stop_pre_postList.group(2))<=int(variation['POS within Protein']) and int(variation['POS within Protein'])<=int(proacc_start_stop_pre_postList.group(3)):
                            ##this peptide is possibly an evidence of the SAP.
                            ## Both vcf location and the protein identification locations (at least for MSGF+)
                            ## are 1 base.But pythonic indexes are 0 based.
                            position=int(variation['POS within Protein'])-int(proacc_start_stop_pre_postList.group(2))
                            if variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position]:
                                ## The peptide supports the SAP/SSAP.
                                ## add this proof to the matrix
                                prevFound=1
                                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Full'})),ignore_index=True)
                    else:
                        if variation['Type']=='ALT' or variation['Type']=='SALT':
                            ## The Alteration event might have partial peptide coverage.
                            altStart=(intvariation['POS within Protein'])
                            altEnd=int(variation['POS within Protein'])+len(variation['ALT'])-1
                            ## Check whether this identified peptide overlaps with this alteration event.
                            ## The overlap can either be full or partial overlap.
                            if int(proacc_start_stop_pre_postList.group(2))<=altStart and altEnd<=int(proacc_start_stop_pre_postList.group(3)):
                                ##this means the alteration event is fully covered by the peptide.
                                ##pythonic indexes are 0 based.
                                position=int(variation['POS within Protein'])-int(proacc_start_stop_pre_postList.group(2))
                                altLength=len(variation['ALT'])
                                if variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position:(position+altLength)]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Full'})),ignore_index=True)
                            elif int(proacc_start_stop_pre_postList.group(2))>=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(2)) and altEnd<=int(proacc_start_stop_pre_postList.group(3)):
                                ## Partial peptide coverage
                                position=???
                                if variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][???]
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Partial'})),ignore_index=True)
                            elif int(proacc_start_stop_pre_postList.group(2))<=altStart and altStart<=int(proacc_start_stop_pre_postList.group(3)) and altEnd>=int(proacc_start_stop_pre_postList.group(3)):
                                ## Partial match
                                position=???
                                if variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][???]
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Partial'})),ignore_index=True)
                        elif variation['Type']=='INS':
                            ## Insertion event
                        elif variation['Type']=='DEL':
                            ## Deletion Event. Need to think how this event can be validated.
            else:
                ##ERROR
                print("ERROR: proacc_start_stop_pre_post list should have single entry")
        else:
            ## This suggest multiple Spectra match to a peptide. Peptide sequence might be the same, but the PSM
            ## is not.
            if prevFound==1:
                ##this means prevPeptide was counted as an evidence of the variaton. Hence this PSM should also
                ##be counted.
def groupPeptide(peptideGrouped, ORFId):
    PSMsList=pd.DataFrame();
    for seq, query in peptideGrouped:
        # here query holds all the PSMs for the peptide Sequence. 'proteinacc_start_stop_pre_post_;' column
        # lists corresponding the proteins and its location within the proteins.
        # Accessing the first row of the group is enough as same sequence will have same protein list.
        protList=query.iloc[0]['proteinacc_start_stop_pre_post_;']
        if ORFId in protList:
            ## this means the ORF in question was identified and containes SAP/ALT/INDELs
            ## This peptide/Sequence is an evidence of the ORF in interest. Hence, put these PSMs together.
            ## Essentially this is grouping PSMs according to proteins. As protein file tells us which
            ## peptides belong to the sub-members of a PAG
            ## create a key value structure, where the orfsid is the key
            PSMsList=PSMsList.append(query,ignore_index=True)
    return PSMsList
  
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
    for name, variations in vcfGrouped:
        ## check whether this ORF has been identified and whether this SAP/ALT/INDEL event has a peptide evidence.
        ## the ORF id contained ',' which was replaced by ';' in contigStat.pl as this code produce comma separated file. Later
        ## on IdentifyProeinIsoformSAP.py replaced ';' by '&' for same cause. ORFs with multiple parent
        ## transcripts have '&'. Whereas, the same ORFs in the PSMs list contain ','.
        pattern=re.compile('&')
        nameComma=pattern.sub(',',name)
        ## find peptide evidence of this ORF
        PSMsOfORF=groupPeptide(peptideGrouped, ORFId)
        ##For each entry in the query, which is essentially the variations, check overlap between the variation and these peptide.
        for index, variation in variations:
            ## Match 
            PSMsOfVariations=checkVariationPeptideCoverage(ORFId, PSMsOfORF, variation)
        
    
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
    
PSMFileName="D:\data\Results\Human-Adeno\Identification\PASA\sORF\pasa_assemblyV1+fdr+th+grouping.csv"
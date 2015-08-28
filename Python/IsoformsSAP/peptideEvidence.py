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
    print("variation id:"+str(variation['ID']))
    for i in range(0,len(PSMsOfORF)):
        
        protAcc=PSMsOfORF.iloc[i]['proteinacc_start_stop_pre_post_;']
        print("prot accession:"+protAcc)
        pepSeq=PSMsOfORF.iloc[i]['Sequence']
        if prevPeptide!=pepSeq:
            ## Check whether this PSM covers the location of the variation.
            ## If yes, compare the peptide sequence and the variation sequence to make sure the peptide does
            ## support the variation.
            ## Split the protId, get start and end and compare with the loc. Check the sequence, compare with altSeq.
            prevPeptide=pepSeq;
            ## prevFound flag is reset here. Because prevPeptide is not same as current one, therefore the PSM has not been
            ## added to the PSMRvidence yet.
            prevFound=0
            protList=protAcc.split(';')
            proacc_start_stop_pre_post=[s for s in protList if ORFId in s] #this should always have one element.
            if len(proacc_start_stop_pre_post)==1:
                ##this should always be the case.
                startEndPrePost=proacc_start_stop_pre_post[0].replace(ORFId,'')
                print("start end pre post:"+startEndPrePost)
                pattern=re.compile('_(\d+)_(\d+)_([A-Z]|-)_([A-Z]|-)')
                match=pattern.finditer(startEndPrePost)
                count=0
                for m in match:
                    ##this loop should run only once.
                    proacc_start_stop_pre_postList=m
                    if count>0:
                        ##ERROR
                        print("ERROR: Match object should have single value."+protAcc)
                        print("ORFId:"+ORFId)
                        print("proacc_start_stop_pre_postList:"+proacc_start_stop_pre_postList.groups())
                        count=count+1;
                        break;
                    count=count+1;
                ## this is what should happen
                #print("Count:"+str(count))
                if count==1:
                    ## Check whether location of the SAP/ALT/INDEL is covered by this peptide.
                    ## For SAPs its straight forward. 
                    ## It might happen that for ALT/INDELs the peptide covers the event partially.
                    #print("Count:"+str(count))
                    if variation['Type']=='SAP' or variation['Type']=='SSAP':
                        #print("Var type S/SAP")
                        print("proacc start-end:"+str(proacc_start_stop_pre_postList.group(1))+"-"+str(proacc_start_stop_pre_postList.group(2)))
                        print("var loc="+variation['QPOS'])
                        if int(proacc_start_stop_pre_postList.group(1))<=int(variation['QPOS']) and int(variation['QPOS'])<=int(proacc_start_stop_pre_postList.group(2)):
                            ##this peptide is possibly an evidence of the SAP.
                            ## Both vcf location and the protein identification locations (at least for MSGF+)
                            ## are 1 base.But pythonic indexes are 0 based.
                            position=int(variation['QPOS'])-int(proacc_start_stop_pre_postList.group(1))
                            print("1. Overlap")
                            if variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position]:
                                ## The peptide supports the SAP/SSAP.
                                ## add this proof to the matrix
                                prevFound=1
                                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Full'})),ignore_index=True)
                                print("1.pep found:")
                                print(PSMEvidence.to_csv())
                            else:
                                print(str(variation['ID'])+" SAP/SSAP event event has peptide covering the position, but AA is not same. peptide :"+PSMsOfORF.iloc[i]['Sequence'])
                        else:
                            print(str(variation['ID'])+": This SAP/SSAP is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
                    else:
                        if variation['Type']=='ALT' or variation['Type']=='SALT':
                            ## The Alteration event might have partial peptide coverage.
                            altStart=int(variation['QPOS'])
                            #print("ORFId:"+ORFId)
                            #print("variation:"+variation.to_csv(None))
                            #print("alt:"+str(variation['ALT']))
                            altEnd=int(variation['QPOS'])+len(variation['ALT'])-1
                            ## Check whether this identified peptide overlaps with this alteration event.
                            ## The overlap can either be full or partial overlap.
                            if int(proacc_start_stop_pre_postList.group(1))<=altStart and altEnd<=int(proacc_start_stop_pre_postList.group(2)):
                                ##this means the alteration event is fully covered by the peptide.
                                ##pythonic indexes are 0 based.
                                position=int(variation['QPOS'])-int(proacc_start_stop_pre_postList.group(1))
                                altLength=len(variation['ALT'])
                                print("2. Overlap")
                                if variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position:(position+altLength)]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Full'})),ignore_index=True)
                                    print("2.pep found")
                                else:
                                    print(str(variation['ID'])+" ALT/SALT event event has peptide covering the position, but AA is not same. peptide :"+PSMsOfORF.iloc[i]['Sequence'])
                            elif int(proacc_start_stop_pre_postList.group(1))>=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(1)) and altEnd<=int(proacc_start_stop_pre_postList.group(2)):
                                ## Partial peptide coverage
                                position=0
                                positionEnd=altEnd-int(proacc_start_stop_pre_postList.group(1))+1
                                altIdx=altStart-int(proacc_start_stop_pre_postList.group(1))
                                print("3. Overlap")
                                if variation['ALT'][altIdx:]==PSMsOfORF.iloc[i]['Sequence'][position:positionEnd]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Partial'})),ignore_index=True)
                                    print("3.pep found")
                                else:
                                    print(str(variation['ID'])+" ALT/SALT event event has peptide covering the position, but AA is not same. peptide :"+PSMsOfORF.iloc[i]['Sequence'])
                            elif int(proacc_start_stop_pre_postList.group(1))<=altStart and altStart<=int(proacc_start_stop_pre_postList.group(2)) and altEnd>=int(proacc_start_stop_pre_postList.group(2)):
                                ## Partial match
                                position=altStart-int(proacc_start_stop_pre_postList.group(1))
                                altEndIdx=int(proacc_start_stop_pre_postList.group(2))-altStart
                                print("4. Overlap")
                                if variation['ALT'][0:altEndIdx]==PSMsOfORF.iloc[i]['Sequence'][position:]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Partial'})),ignore_index=True)
                                    print("4.pep found")
                                else:
                                    print(str(variation['ID'])+":This ALT/SALT event has peptide covering the position, but AA is not same. Peptide :"+PSMsOfORF.iloc[i]['Sequence'])
                            else:
                                print(str(variation['ID'])+": This variation is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
                        elif variation['Type']=='INS':
                            ## Insertion event. This might have partial or full peptide coverage as the alteration events.
                            altStart=int(variation['QPOS'])
                            altEnd=int(variation['QPOS'])+len(variation['ALT'])-1
                            ## Check whether this identified peptide overlaps with this insertion event.
                            ## The overlap can either be full or partial overlap.
                            if int(proacc_start_stop_pre_postList.group(1))<=altStart and altEnd<=int(proacc_start_stop_pre_postList.group(2)):
                                ##this means the insertion event is fully covered by the peptide.
                                ##pythonic indexes are 0 based.
                                position=int(variation['QPOS'])-int(proacc_start_stop_pre_postList.group(1))
                                altLength=len(variation['ALT'])
                                print("5. Overlap")
                                if variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position:(position+altLength)]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Full'})),ignore_index=True)
                                    print("5.pep found")
                                else:
                                    print(str(variation['ID'])+":This INS event has peptide covering the position, but AA is not same. Peptide :"+PSMsOfORF.iloc[i]['Sequence'])
                            elif int(proacc_start_stop_pre_postList.group(1))>=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(1)) and altEnd<=int(proacc_start_stop_pre_postList.group(2)):
                                ## Partial peptide coverage
                                position=0
                                positionEnd=altEnd-int(proacc_start_stop_pre_postList.group(1))+1
                                ## part of the variation['ALT'] will be missing
                                altIdx=altStart-int(proacc_start_stop_pre_postList.group(1))
                                print("6. Overlap")
                                if variation['ALT'][altIdx:]==PSMsOfORF.iloc[i]['Sequence'][position:positionEnd]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Partial'})),ignore_index=True)
                                    print("6.pep found")
                                else:
                                    print(str(variation['ID'])+": This INS event has peptide covering the position, but AA is not same. Peptide :"+PSMsOfORF.iloc[i]['Sequence'])
                            elif int(proacc_start_stop_pre_postList.group(1))<=altStart and altStart<=int(proacc_start_stop_pre_postList.group(2)) and altEnd>=int(proacc_start_stop_pre_postList.group(2)):
                                ## Partial match
                                position=altStart-int(proacc_start_stop_pre_postList.group(1))
                                altEndIdx=int(proacc_start_stop_pre_postList.group(2))-altStart+1 ##1 added because this is second index. In python, right hand limit is exclusive
                                print("7. Overlap")
                                if variation['ALT'][0:altEndIdx]==PSMsOfORF.iloc[i]['Sequence'][position:]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Partial'})),ignore_index=True)
                                    print("7.pep found")
                                else:
                                    print(str(variation['ID'])+": This INS event has peptide covering the position, but AA is not same. Peptide :"+PSMsOfORF.iloc[i]['Sequence'])
                            else:
                                print(str(variation['ID'])+": This variation is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
                        elif variation['Type']=='DEL':
                            ## Deletion Event. To proof that the deletion is true, I need to check that the peptide sequence does not support the
                            ## ref sequence.
                            altStart=int(variation['QPOS'])
                            altEnd=int(variation['QPOS'])+len(variation['REF'])-1
                            ## Check whether this identified peptide overlaps with this deletion event.
                            ## The overlap can either be full or partial overlap.
                            if int(proacc_start_stop_pre_postList.group(1))<=altStart and altEnd<=int(proacc_start_stop_pre_postList.group(2)):
                                ##this means the deletion event is fully covered by the peptide.
                                ##pythonic indexes are 0 based.
                                position=int(variation['QPOS'])-int(proacc_start_stop_pre_postList.group(1))
                                altLength=len(variation['REF'])
                                print("8. Overlap")
                                ## if the ref sequenece is there in the peptide or the only the deleted AAs are there, then the deletion event is not supported.
                                ## The OR part might be true due to AAs repeatation, but we dont have a way to disproof it without further checking. 
                                if variation['REF']==PSMsOfORF.iloc[i]['Sequence'][position:(position+altLength)] or variation['REF'][1:]==PSMsOfORF.iloc[i]['Sequence'][(position+1):(position+altLength)]:
                                    print(str(variation['ID'])+" DEL event event has peptide covering the position. The peptide sequence supports the reference sequence :"+PSMsOfORF.iloc[i]['Sequence'])
                                elif variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position:(position+1)]:
                                    ## position plus one because alt will only have one AA from deletion start point-1 position to represent the deletion.
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'Full'})),ignore_index=True)
                                    print("8.pep found")
                                else:
                                    print("ERROR: peptide neither support ALT nor REF. variation:"+variation['ID']+", peptide:"+PSMsOfORF.iloc[i]['Sequence'])
                            elif int(proacc_start_stop_pre_postList.group(1))>=altStart and altEnd>=int(proacc_start_stop_pre_postList.group(1)) and altEnd<=int(proacc_start_stop_pre_postList.group(2)):
                                ## Partial peptide coverage
                                position=0
                                positionEnd=altEnd-int(proacc_start_stop_pre_postList.group(1))+1
                                altIdx=altStart-int(proacc_start_stop_pre_postList.group(1))
                                print("9. Overlap")
                                if variation['REF'][altIdx:]==PSMsOfORF.iloc[i]['Sequence'][position:positionEnd] or variation['REF'][(altIdx+1):]==PSMsOfORF.iloc[i]['Sequence'][(position+1):positionEnd]:
                                    print(str(variation['ID'])+" DEL event event has peptide covering the position, but the peptide is partially supporting the reference sequence:"+PSMsOfORF.iloc[i]['Sequence'])
                                elif variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'FULL'})),ignore_index=True)
                                    print("9.pep found")
                                else:
                                    print("ERROR: peptide neither support ALT nor REF. variation:"+str(variation['ID'])+", peptide:"+PSMsOfORF.iloc[i]['Sequence'])
                            elif int(proacc_start_stop_pre_postList.group(1))<=altStart and altStart<=int(proacc_start_stop_pre_postList.group(2)) and altEnd>=int(proacc_start_stop_pre_postList.group(2)):
                                ## Partial match
                                position=altStart-int(proacc_start_stop_pre_postList.group(1))
                                altEndIdx=int(proacc_start_stop_pre_postList.group(2))-altStart+1
                                print("10. Overlap")
                                if variation['REF'][0:altEndIdx]==PSMsOfORF.iloc[i]['Sequence'][position:] or variation['REF'][1:altEndIdx]==PSMsOfORF.iloc[i]['Sequence'][(position+1):]:
                                    print(str(variation['ID'])+" DEL event event has peptide covering the position, but the peptide is partially supporting the reference sequence:"+PSMsOfORF.iloc[i]['Sequence'])
                                elif variation['ALT']==PSMsOfORF.iloc[i]['Sequence'][position]:
                                    prevFound=1
                                    PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':'FULL'})),ignore_index=True)
                                    print("10.pep found")
                                else:
                                    print("ERROR: peptide neither support ALT nor REF. variation:"+str(variation['ID'])+", peptide:"+PSMsOfORF.iloc[i]['Sequence'])
                            else:
                                print(str(variation['ID'])+": This variation is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
                elif count==0:
                    print("Count is 0, i.e. protacc_start_end_prev_next does not match:"+proacc_start_stop_pre_post[0])
                    print("ORFId:"+ORFId)
                else:
                    print("ERROR: count is not 1:"+protAcc)
                    print("ORFId:"+ORFId)
                    print("proacc_start_stop_pre_postList:"+proacc_start_stop_pre_postList.groups())
            else:
                ##ERROR
                print(proacc_start_stop_pre_post+" ERROR: proacc_start_stop_pre_post list should have single entry")
        else:
            ## This suggest multiple Spectra match to a peptide. Peptide sequence might be the same, but the PSM
            ## is not.
            if prevFound==1:
                print("11. Overlap")
                print("10.pep found")
                print("Length of PSMEvidence"+str(len(PSMEvidence)))
                ##this means prevPeptide was counted as an evidence of the variaton. Hence this PSM should also
                ##be counted.
                print(PSMsOfORF.iloc[i].to_csv(None))
                PSMEvidence=PSMEvidence.append(PSMsOfORF.iloc[i].append(pd.Series({'Evidence':PSMEvidence.iloc[len(PSMEvidence)-1]['Evidence']})),ignore_index=True)
            #else:
            #    print(str(variation['ID'])+": This variation is not supported by "+PSMsOfORF.iloc[i]['Sequence'])
    print("PSMEvd:"+PSMEvidence.to_csv())
    return PSMEvidence

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

def addVariationInfoToPSM(PSMsOfVariations, variation):
    ## Add ORFId and variation information.
    varInf={'ORFID':variation['QueryID'],'SubjectID':variation['SubjectID'],'variationID':variation['ID'],'variationType':variation['Type'],'Location':variation['QPOS']}
    lenPSMVar=len(PSMsOfVariations)
    writeDF=pd.concat([PSMsOfVariations,pd.DataFrame([varInf]*lenPSMVar,index=range(lenPSMVar))],axis=1)
    return writeDF

def printPSMsOfVariation(writeDF, newPSMFile, headerFlag):
    ##this is PSM wise view
    if headerFlag==1:
        ## Write column names
        writeDF.to_csv(newPSMFile, index=False)
    else:
        writeDF.to_csv(newPSMFile, header=False, index=False)

def calculateQualityScore(PSMsOfVariation,featureWeights):
    ## For time being this function returns 1, but in future I need to figure out a formula to represent the confidence of peptide evidences supporting a variation.
    return 1

def addPSMInfoToVariation(variation, PSMsOfVariation):
    ##Add PSM count, Peptide count, possibly an average score indicating PSM/peptide quality.
    ##For time being quality score is 1.
    score=calculateQualityScore(PSMsOfVariation,1)
    ## Number of PSMs
    PSMCount=len(PSMsOfVariation['PSM_ID'].unique())
    ## Number of peptide sequences
    peptideCount=len(PSMsOfVariation['Sequence'].unique())
    
    ## Count of peptide sequences that uniquely map to this ORF/Protein
    unqPeptideCount=len(PSMsOfVariation[(~PSMsOfVariation['proteinacc_start_stop_pre_post_;'].str.contains(';'))]['Sequence'].unique()) ##~ is doing the job of negation, i.e. does not contain. This means the peptide is not shared with other proteins
    ## Unique peptides
    print("Unique Peptide:")
    print(PSMsOfVariation[(~PSMsOfVariation['proteinacc_start_stop_pre_post_;'].str.contains(';'))]['Sequence'].unique())
    unqPeptides=",".join(PSMsOfVariation[(~PSMsOfVariation['proteinacc_start_stop_pre_post_;'].str.contains(';'))]['Sequence'].unique())
    if not unqPeptides:
        unqPeptides="-"
    PSMInfo={'PSMCount':PSMCount,'PeptideCount':peptideCount,'UniquePeptideCount':unqPeptideCount,'UniquePeptide':unqPeptides,'ConfidenceScore':score}
    writeVar=variation.append(pd.Series(PSMInfo))
    print(writeVar.to_csv(None))
    return writeVar

def printValidatedVariations(variation, newVcfFile, headerFlag):
    ##This is another vcf
    if headerFlag==1:
        newVcfFile.write("#Chr\tPOS within Protein\tID\tREF\tALT\tINFO\n")
    ##SubjectID','QueryID','Alignment','Type','QPOS'
    newVcfFile.write(str(variation['#Chr'])+"\t"+str(variation['POS within Protein'])+"\t"+str(variation['ID'])+"\t"+variation['REF']+"\t"+variation['ALT']+"\t"+"SubjectId="+variation['SubjectID']+";QueryId="+variation['QueryID'].replace(',','&')+";Alignment=["+variation['Alignment']+"];Type:"+variation['Type']+";QPOS:"+str(variation['QPOS'])+";PeptideCount:"+str(variation['PeptideCount'])+";UniquePeptideCount:"+str(variation['UniquePeptideCount'])+";Peptides:"+variation['UniquePeptide']+";Score:"+str(variation['ConfidenceScore'])+"\n")

def findPeptideEvidence(vcf, PSMs, newVcfFile, newPSMFile):
    ## vcf file fields: Chr, POS within Protein,ID, REF, ALT, INFO(SubjectId=P09417-2;QueryId=Dataset_A_asmbl_41426_ORF20_Frame_3_84-446;Alignment=[QueryLength=121:QueryStart=1:QueryEnd=116:SubjectLength=213:SubjectStart=1:SubjectEnd=116];Type:SSAP;QPOS:116)
    ## split the INFO column and add them to main data frame
    info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    print(vcf.QueryID[0:5])
    vcf.Alignment=vcf.Alignment.str.replace('Alignment=\[','')
    vcf.Alignment=vcf.Alignment.str.replace('\]','')
    vcf.Type=vcf.Type.str.replace('Type:','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS:','')
    ##groups vcf entries by ORF/Query ids.
    vcfGrouped=vcf.groupby('QueryID')
    
    ##in the same way group the PSMs according to the prptide sequence.
    peptideGrouped=PSMs.groupby('Sequence')
    # each of these group represents all the SAPs/ALTs/INDELs for a ORF/Query
    headerFlag=1
    for name, variations in vcfGrouped:
        ## check whether this ORF has been identified and whether this SAP/ALT/INDEL event has a peptide evidence.
        ## the ORF id contained ',' which was replaced by ';' in contigStat.pl as this code produce comma separated file. Later
        ## on IdentifyProeinIsoformSAP.py replaced ';' by '&' for same cause. ORFs with multiple parent
        ## transcripts have '&'. Whereas, the same ORFs in the PSMs list contain ','.
        pattern=re.compile('&')
        ORFId=pattern.sub(',',name)
        ## find peptide evidence of this ORF
        print("ORFId:"+ORFId)
        
        PSMsOfORF=groupPeptide(peptideGrouped, ORFId)
        if len(PSMsOfORF)>0:    
            ##For each entry in the query, which is essentially the variations, check overlap between the variation and these peptide.
            for i in range(len(variations)):
                variation=variations.iloc[i]
                ## Match
                PSMsOfVariations=checkVariationPeptideCoverage(ORFId, PSMsOfORF, variation)
                #print(PSMsOfVariations.head(1))
                if len(PSMsOfVariations)>0:
                    print("Variation "+str(variation['ID'])+" has peptide evidence")
                    PSMsOfVariations=addVariationInfoToPSM(PSMsOfVariations, variation)
                    variation=addPSMInfoToVariation(variation, PSMsOfVariations)
                    printPSMsOfVariation(PSMsOfVariations, newPSMFile,headerFlag)
                    printValidatedVariations(variation, newVcfFile,headerFlag)
                    headerFlag=0
                else:
                    print(str(variation['ID'])+" does not have peptide evidence")
        else:
            print(ORFId+" was not identified")
        

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def main(PSMFileName, vcfFileName, newVcfFileName, newPSMFileName):
    ##read peptide identification file
    PSMs=readFile(PSMFileName, ',')
    ##read vcf file
    vcf=readFile(vcfFileName,'\t')
    ##for each entry in the vcf try find peptides that overlaps the vcf entry location. For deletion event, it might
    ##be little tricky. Not finding any peptide for deletion event is a good sign but does not confirms the deletion.
    with open(newVcfFileName, 'w') as newVcfFile, open(newPSMFileName, 'w') as newPSMFile:
        findPeptideEvidence(vcf, PSMs, newVcfFile, newPSMFile)


PSMFileName="D:/data/Results/Human-Adeno/Identification/PASA/sORF/pasa_assemblyV1+fdr+th+grouping.csv"
vcfFileName="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_with_Location_VariationV7.vcf"

## New Files
newPSMFileName="D:/data/Results/Human-Adeno/Identification/PASA/sORF/pasa_assemblyV1+fdr+th+groupingVariationEvidence.csv"
newVcfFileName="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_with_Location_VariationV7PeptideEvidence.vcf"
'''
PSMFileName="pepTest.csv"
vcfFileName="vcfTest.vcf"

## New Files
newPSMFileName="pepVariationEvidence.csv"
newVcfFileName="vcfPeptideEvidence.vcf"
'''
main(PSMFileName, vcfFileName, newVcfFileName, newPSMFileName)
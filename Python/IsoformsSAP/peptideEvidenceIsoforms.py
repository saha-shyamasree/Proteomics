## this code finds peptide evidence of isoforms.

import pandas as pd
import argparse
from AminoAcidVariation import AminoAcidVariation

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj;

def filterPeptide(peptideObj, protIds, pepColumn):
    peptideObj['Is decoy']=peptideObj['Is decoy'].astype(str)
    #print(peptideObj['Is decoy'])
    protStr="|".join(protIds)
    ##removing decoy hits
    peptideObj=peptideObj[~peptideObj['Is decoy'].str.contains("TRUE",case=False)]
    ##removing PSMs only mapping to Contaminents
    peptideObj=peptideObj[~peptideObj[pepColumn].str.contains("^(CONT.*;?)+$")]
    ##removing PSMs only mapping to Decoy
    peptideObj=peptideObj[~peptideObj[pepColumn].str.contains("^(XXX.*;?)+$")]
    filteredPeptideObj = peptideObj[peptideObj[pepColumn].str.contains(protStr)]
    return filteredPeptideObj

def pepThresholding(prots, pepTh, protRevStr, protContStr):
    print("1. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protRevStr)]
    print("2. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains(protContStr)]
    print("3. prots dim:"+str(prots.shape))
    prots=prots[prots['distinct peptide sequences']>pepTh]
    print("4. prots dim:"+str(prots.shape))
    return prots

def filterBlast(blastFile, queryList):
    blast=readFile(blastFile,",")
    blast['query_name'].isin(queryList)
    identIdx=blast['query_name'].isin(queryList)
    idx=identIdx[identIdx==True].index.tolist()
    blastIdentified = blast.iloc[idx,:]
    return blastIdentified

def createIsoObject(pos,qpos,isoClass, varCount, ref, alt, iso, count):
    varCount=varCount+1
    alignInfo={'queryStart':iso['q_st'],'queryEnd':iso['q_end'],'qLen':iso['query_length'],'subjectStart':iso['s_st'],'subjectEnd':iso['s_end'],'sLen':iso['hit_length'],'qSerialNo':count}
    iso=AminoAcidVariation(iso['hit_def'],iso['query_name'],pos,qpos,ref,alt,isoClass,iso['Location'],varCount,alignInfo)
    return (varCount,count, iso)

def findORFType(header):
    orfType=''
    if 'type:complete' in header:
        orfType='_complete'
    elif 'type:3prime_partial' in header:
        orfType='_C-terminus_partial'
    elif 'type:5prime_partial' in header:
        orfType='_N-terminus_partial'
    elif 'type:internal' in header:
        orfType='_internal'
    else:
        print("Error: ORF type was not found. This should not happen if the ORFs were predicted using transdecoder.")
    return orfType

def IsoformVariation(identifiedBlast,count):
    ##we dont want to store the sequence, hence alt and ref are both '.'
    alt="."
    ref="."
    ##positions need redoing.
    isoList=[]
    for i,iso in identifiedBlast.iterrows():
        varCount=0
        count=count+1
        orfType=findORFType(iso['query_name'])
        if iso['q_st']==1:
            if iso['s_st']==1:
                if iso['q_end']==iso['query_length']:
                    if iso['s_end']==iso['hit_length']:
                        print("This is a case of known or known with variation protein. This ORF should not be in the list of isoforms.")
                    elif iso['s_end']<iso['hit_length']:
                        #The protein is partially mapped to the ORF. Because the q_end==query_length, we can say that the protein is longer than the ORF.
                        #The ORF is shorter than the protein at the 3' end.
                        isoClass='DEL_3prime'+orfType   
                        varCount, count, isoObj = createIsoObject(iso['s_end'],iso['q_end'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    if iso['s_end']==iso['hit_length']:
                        ##The ORF is larger than the protein, because the full length of the protein has a map  to the ORF.
                        isoClass='INS_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        ##Both the ORF and the protein has sequence at the 3' end that does not map to each other. Might be result of alternative
                        ## end exon.
                        isoClass='ALT_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print('queary_end can not be bigger than the query length')
            elif iso['s_st']>1:
                ##protein pos is sent as 1 because the AminoAcid variation print function deals with iso['s_st']>1 ans add the value as prefix.
                if iso['q_end']==iso['query_length']:
                    if iso['s_end']==iso['hit_length']:
                        ##protein's 5' side has extra amino acids. i.e. this is a deletion event at 5'
                        isoClass='DEL_5prime'+orfType
                        varCount, count, isoObj=createIsoObject(1, iso['q_st'], isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        ##part of the protein match with the ORF. But the whole ORF is mapping to the protein. So two events.
                        isoClass='DEL_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(1,iso['q_st'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='DEL_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1, iso['q_end'], isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    #q_end<quety_length, i.e. the whole ORF did not match.
                    if iso['s_end']==iso['hit_length']:
                        ##this is an event of INS at the 3' end. And also an event of deletion at the 5' end. because s_st>1
                        isoClass='DEL_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(1,iso['q_st'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='INS_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        ##Both the ORF and the protein has sequence at the 3' end that does not map to each other. Might be result of alternative
                        ## end exon.
                        isoClass='DEL_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(1,iso['q_st'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        orfType=findORFType(iso['query_name'])
                        isoClass='ALT_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end'],isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print('queary_end can not be bigger than the query length')
        elif iso['q_st']>1:
            if iso['s_st']==1:
                ##Insertion
                if iso['q_end']==iso['query_length']:
                    #3' side ORF complete match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        isoClass='INS_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3'del
                        isoClass='INS_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='DEL_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    #3' side ORF short match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        #as q_st>1 and q_end<q_length, we have insertion at both the end.
                        isoClass='INS_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='INS_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3' ALT
                        isoClass='INS_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_st'],1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='ALT_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end'],iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print("Error: query_end can not be bigger than the query length")
            elif iso['s_st']>1:
                ##check Updated code. Here both query and subject 5 prime entries shouls have position value 1.
                ##Alteration at 5 prime.
                if iso['q_end']==iso['query_length']:
                    #3' side ORF complete match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        isoClass='ALT_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3'del
                        isoClass='ALT_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='DEL_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                elif iso['q_end']<iso['query_length']:
                    #3' side ORF short match
                    if iso['s_end']==iso['hit_length']:
                        #3' protein full match
                        #as q_st>1 and q_end<q_length, we have insertion at both the end.
                        isoClass='ALT_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='INS_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    elif iso['s_end']<iso['hit_length']:
                        #3' protein short match
                        #so event of 5'INS and 3' ALT
                        isoClass='ALT_5prime'+orfType   
                        varCount, count, isoObj=createIsoObject(1,1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                        isoClass='ALT_3prime'+orfType   
                        varCount, count, isoObj=createIsoObject(iso['s_end']-iso['s_st']+1,iso['q_end']-iso['q_st']+1,isoClass, varCount, ref, alt, iso, count)
                        isoList.append(isoObj)
                    else:
                        print("Error: subject_end can not be bigger than the subject length")
                else:
                    print("Error: query_end can not be bigger than the query length")
                
    return (count, isoList)


def checkIsoformPeptideCoverage(ORFId, PSMsOfORF, variation):
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
                ## 'start' and 'stop' are both inclusive
                startEndPrePost=proacc_start_stop_pre_post[0].replace(ORFId,'')
                print("start end pre post:"+startEndPrePost)
                pattern=re.compile('_(\d+)_(\d+)_([A-Z]|-)_([A-Z]|-)')
                match=pattern.finditer(startEndPrePost)
                count=0

def groupPeptide(psms, ORFId, pepColumn, rev):
    ##not considering CONT because CONT entries are not going to have ORFId.
    regex="(?<!"+rev+")"+re.escape(ORFId)
    filteredPSMs= psms[psms[pepColumn].str.contains(regex)]
    return filteredPSMs

def groupVariation(variationList, ORFId):
    ##not considering CONT because CONT entries are not going to have ORFId.
    regex="(?<!"+rev+")"+re.escape(ORFId)
    filteredVar= psms[psms[pepColumn].str.contains(regex)]
    return filteredPSMs

##isoList is a list of AminoAcidVariation class.
def peptidePerIsoform(peptides, variationList, isoIds, isoList, pepColumn, rev):
    print("Peptide Count:")
    ##Known proteins
    psms=filterPeptide(peptides, isoIds, pepColumn)
    #psms.to_csv(,index=False)
    
    print("Peptides:"+str(len(psms['Sequence'].unique())))
    for i in range(0,len(isoIds)):
        PSMsOfORF=groupPeptide(psms, isoIds[i], pepColumn, rev)
        variations=groupVariation(variationList, isoIds[i])
        checkIsoformPeptideCoverage(isoIds[i], PSMsOfORF, variations)

 
    
protRevStr="XXX_"
protContStr="CONT_"
pepTh=1
pepColumn='proteinacc_start_stop_pre_post_;'
parser = argparse.ArgumentParser(description='This script read output of UniProteinLocation.py and identify variations')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="full path of blast csv file", metavar="PATH")
parser.add_argument("-i", "--isoform", nargs=1, required=True, help="full path to the isoform file", metavar="PATH")
#parser.add_argument("-j", "--isovar", nargs=1, required=True, help="full path to the isoforms with variation file", metavar="PATH")
#parser.add_argument("-p", "--pep", nargs=1, required=True, help="full path to the peptide identification file", metavar="PATH")
#parser.add_argument("-o", "--output", nargs=1, required=True, help="full path to the output csv file", metavar="PATH")
args = parser.parse_args()

##Read iso file
isoforms=readFile(args.isoform[0],',')
### SEparates blast results for the isoform list
identifiedIsoBlast=filterBlast(args.blast[0], isoforms['ORF Id'])

count=0
count, isoList=IsoformVariation(identifiedIsoBlast,count)

for i in range(0,len(isoList)):
    print(isoList[i].toString())
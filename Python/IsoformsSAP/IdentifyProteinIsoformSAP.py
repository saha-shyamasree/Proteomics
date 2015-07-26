###This code identifies SAPs/Alteration/INDELs in an ORF BLAST map to an Uniprot protein.
###While identifying the SAPs/Alt/INDELs, it produces a list of exactly match between ORFs
###and uniprot proteins that is known proteins, a list of known proteins with SAPs/Alt/INDELs
###and a list of possible protein isoforms. These lists can then easily be read into R to check how
###many of these ORFs are identified from MS experient.
import os
import csv
import re
from AminoAcidVariation import AminoAcidVariation

def checkMatch(data, value):
    if data==value:
        return 1;
    else:
        return 0;

def floatMatchCheck(data, value):
    if data==value:
        return 1;
    else:
        return 0;

def floatAboveThreshold(data, value):
    if data>=value:
        return 1;
    else:
        return 0;

def evalueCheck(data, threshold):
    if data<=threshold:
        return 1;
    else:
        return 0;

def lengthRatioCheck(qLen, sLen):
    if float(qLen)/sLen==1:
        return 1;
    else:
        return 0;

def synonimousSAPs(seq1,seq2,seq3,delCount,varCount,aaVariationList,pos,sid,qid,chrm, alignInfo): #seq1 is query, seq2 is ref, seq3 is match from blast.
    match=re.finditer('\++',seq3);
    for m in match:
        if len(m.group(0))==1:
            aaVariationList.append(AminoAcidVariation(sid, qid, pos+m.start(0)+1-delCount, seq2[m.start(0):m.end(0)], seq1[m.start(0):m.end(0)], "SSAP", chrm, varCount,alignInfo))
        else:
            aaVariationList.append(AminoAcidVariation(sid, qid, pos+m.start(0)+1-delCount, seq2[m.start(0):m.end(0)], seq1[m.start(0):m.end(0)], "SALT", chrm, varCount,alignInfo))
        varCount=varCount+1;
    return {'varList':aaVariationList,'varCount':varCount}
    

def printVariationList(varList):
    for l in varList:
        l.printVar();
'''
def findInsertion(d, prevD, ref, alt, aaVariationList, insertionCount, sId, qId, pos, chrm, varCount, alignInfo):
    if prevD==None:
        pos1=pos+1;
        st=0
    else:
        pos1=pos+1+prevD.end(0)-insertionCount;
        st=prevD.end(0);
    altPart1=alt[st:d.start(0)];
    refPart1=ref[st:d.start(0)];
    if "-" not in refPart1: #it means this is an event of SAP/Alteration.
        if len(refPart1)==1:
            aaVariationList.append(AminoAcidVariation(sId,qId,pos1,refPart1,altPart1,"SAP",chrm,varCount,alignInfo))
            varCount=varCount+1;
        else:
            aaVariationList.append(AminoAcidVariation(sId,qId,pos1,refPart1,altPart1,"ALT",chrm,varCount,alignInfo))
            varCount=varCount+1;
        ##adding the deletion
        refPart2=ref[(d.start(0)-1):d.end(0)]
        altPart2=alt[d.start(0)-1]
        aaVariationList.append(AminoAcidVariation(sId,qId,pos+d.end(0)-insertionCount,refPart2,altPart2,"DEL",chrm,varCount,alignInfo))
        varCount=varCount+1;
    else: #this means insertion event has occured.
        ins=re.finditer('\-+',refPart1)
        prevI=None;
        for i in ins:
            if i.start(0)!=0:
                if prevD==None:#prevd is checked cause, if prevd and prevI are both none then default value of pos1=pos+1 is true.
                    if prevI!=None:
                        pos1=pos+1+prevI.end(0)-insertionCount
                        st=prevI.end(0);
                    else:
                        #default values from above;
                        continue;
                else:
                    if prevI!=None:
                        pos1=pos+1+prevD.end(0)+prevI.end(0)-insertionCount;
                        st=prevD.end(0)+prevI.end(0);
                    else:
                        pos1=pos+1+prevD.end(0);
                        st=prevD.end(0)+prevI.end(0);
                if prevI==None and prevD!=None: 
                    pos1=pos+1+prevD.end(0)-insertionCount
                    st=pos+prevD.end(0)
 

                aPart1=altPart1[st:i.start(0)]
                rPart1=refPart1[st:i.start(0)]
                if len(rPart1)==1:
                    aaVariationList.append(AminoAcidVariation(sId,qId,pos1,rPart1,aPart1,"SAP",chrm,varCount,alignInfo))
                    varCount=varCount+1;
                else:
                    aaVariationList.append(AminoAcidVariation(sId,qId,pos1,rPart1,aPart1,"ALT",chrm,varCount,alignInfo))
                    varCount=varCount+1;
                #adding the insertion.
                rPart2=ref[(i.start(0)-1)]
                aPart2=alt[i.start(0)-1:i.end(0)]
                aaVariationList.append(AminoAcidVariation(sId,qId,pos+prevD.end(0)+i.start(0)-insertionCount,rPart2,aPart2,"INS",chrm,varCount,alignInfo))
                varCount=varCount+1;
                insertionCount=insertionCount+len(i.group(0));
                prevI=i
            else:#insertion from the first position.
                #need review
                rPart=refPart1[0]
                aPart=altPart1[0:i.end(0)]
                aaVariationList.append(AminoAcidVariation(sId,qId,pos+prevD.end(0)+i.start(0)-insertionCount,rPart,aPart,"INS",chrm,varCount,alignInfo))
                varCount=varCount+1;
                insertionCount=insertionCount+len(i.group(0));
                prevI=i
        if prevI!=None:# this should be true as refPart1 contains '-'
            if prevI.end(0)<len(alt): # this ensures that insertion event did not occur at the boundary of alt. If it did, the variation has already been taken care of in above loop.
                altPart=alt[prevI.end(0):];
                refPart=ref[prevI.end(0):];
                if len(refPart)==1:
                    aaVariationList.append(AminoAcidVariation(sId,qId,pos1,refPart,altPart,"SAP",chrm,varCount,alignInfo))
                    varCount=varCount+1;
                else:
                    aaVariationList.append(AminoAcidVariation(sId,qId,pos1,refPart,altPart,"ALT",chrm,varCount,alignInfo))
                    varCount=varCount+1;
    return {'varList':aaVariationList, 'varCount':varCount, 'instCount':insertionCount}

def indelSeparation(alt,ref,varCount,aaVariationList,pos,sid,qid,chrm,alignInfo):
    #example alt=ABD---DV-TG and ref=BDATFDT-GAT
    dels=re.finditer('\-+',alt)
    prevD=None;
    insertionCount=0;
    for d in dels:
        if d.start(0)!=0:
            if prevD==None:
                pos1=pos+1;
                st=0
            else:
                pos1=pos+1+prevD.end(0)-insertionCount;
                st=prevD.end(0);
            altPart1=alt[st:d.start(0)];
            refPart1=ref[st:d.start(0)];
            if "-" not in refPart1: #it means this is an event of SAP/Alteration.
                if len(refPart1)==1:
                    aaVariationList.append(AminoAcidVariation(sId,qId,pos1,refPart1,altPart1,"SAP",chrm,varCount,alignInfo))
                    varCount=varCount+1;
                else:
                    aaVariationList.append(AminoAcidVariation(sId,qId,pos1,refPart1,altPart1,"ALT",chrm,varCount,alignInfo))
                    varCount=varCount+1;
                ##adding the deletion
                refPart2=ref[(d.start(0)-1):d.end(0)]
                altPart2=alt[d.start(0)-1]
                aaVariationList.append(AminoAcidVariation(sId,qId,pos+d.end(0)-insertionCount,refPart2,altPart2,"DEL",chrm,varCount,alignInfo))
                varCount=varCount+1;
            else: #this means insertion event has occured.
                ins=re.finditer('\-+',refPart1)
                prevI=None;
                for i in ins:
                    if i.start(0)!=0:
                        if prevI==None and prevD!=None: #prevd is checked cause, if prevd is not none then pos+1 is not true.
                            pos1=pos+1+prevD.end(0)-insertionCount
                            st=pos+prevD.end(0)
                        elif prevI!=None:
                            pos1=pos+1+prevD.end(0)+prevI.end(0)-insertionCount;
                            st=pos+prevD.end(0)+prevI.end(0);
                        aPart1=altPart1[st:i.start(0)]
                        rPart1=refPart1[st:i.start(0)]
                        if len(refPart1)==1:
                            aaVariationList.append(AminoAcidVariation(sId,qId,pos1,rPart1,aPart1,"SAP",chrm,varCount,alignInfo))
                            varCount=varCount+1;
                        else:
                            aaVariationList.append(AminoAcidVariation(sId,qId,pos1,rPart1,aPart1,"ALT",chrm,varCount,alignInfo))
                            varCount=varCount+1;
                        #adding the insertion.
                        rPart2=ref[(i.start(0)-1)]
                        aPart2=alt[i.start(0)-1:i.end(0)]
                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+prevD.end(0)+i.start(0)-insertionCount,rPart2,aPart2,"INS",chrm,varCount,alignInfo))
                        varCount=varCount+1;
                        insertionCount=insertionCount+len(i.group(0));
                        prevI=i
                    else:#insertion from the first position.
                        #need review
                        rPart=refPart1[0]
                        aPart=altPart1[0:i.end(0)]
                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+prevD.end(0)+i.start(0)-insertionCount,rPart,aPart,"INS",chrm,varCount,alignInfo))
                        varCount=varCount+1;
                        insertionCount=insertionCount+len(i.group(0));
                        prevI=i
        else: #deletion from the first position.##review, need ref[-1]
            refPart=ref[0:d.end(0)]
            altPart=alt[0]
            aaVariationList.append(AminoAcidVariation(sId,qId,pos+d.end(0)-insertionCount,refPart2,altPart2,"DEL",chrm,varCount,alignInfo))
            varCount=varCount+1;
        prevD=d;
    #check last bit for insertion, it should be same as from line no 59. can i put that in a function?
   ''' 

def printVariationListtoFile(varList, fileWrt):
    print("writing to file")
    for l in varList:
        fileWrt.write(l.toString()+"\n");

##Genomic/proteomic locations are 1 based, but pythonic indexes are 0 based. Hence, while converting the pythonic location to genomic/proteomic location, it has been adjusted.
def findSAPsAndINDELs(qSeq, sSeq, mSeq, sId, qId, chromosome, alignInfo):
    aaVariationList=[]
    mixCount=0;
    if len(qSeq)==len(sSeq) and len(sSeq)==len(mSeq):
        ##this should always be true.
        match=re.finditer('\s+',mSeq);
        prevM=None;
        varCount=1;
        deletionCount=0;
        seqSt=0;
        for m in match:
            ##+ represents alternative amino acid. hence, that is part of alteration, not deletion. Whereas space may mean deletion or insertion or alteration from the query.
            if prevM==None:
                qSeqPart=qSeq[0:m.start(0)];
                sSeqPart=sSeq[0:m.start(0)];
                mSeqPart=mSeq[0:m.start(0)];
                pos=m.start(0); ## 1 is not added as we want the position before the begining of the INDEL.
                seqSt=0; # position is made 1 base in the synonimousSAPs function
            else:
                qSeqPart=qSeq[prevM.end(0):m.start(0)];
                sSeqPart=sSeq[prevM.end(0):m.start(0)];
                mSeqPart=mSeq[prevM.end(0):m.start(0)];
                pos=m.start(0)-deletionCount;# 1 is not added as we want the position before the begining of the INDEL.
                seqSt=prevM.end(0);
            ##returns an array of SAPs, i.e. AminoAcidVariation objects and latest count of variations to populate 'ID' for the vcf like file.
            retRes=synonimousSAPs(qSeqPart,sSeqPart,mSeqPart,deletionCount,varCount,aaVariationList,seqSt,sId,qId,chromosome, alignInfo);
            aaVariationList=retRes['varList'];
            varCount=retRes['varCount'];
            span=m.end(0)-m.start(0);
            if sSeq[m.start(0):m.end(0)] == ("-"*span): # this means, its an event of insertion
                deletionCount=deletionCount+span;
                if m.start(0)!=0:
                    ref=sSeq[m.start(0)-1] # -1 to get the index of the amino acid before the deletion location.
                    alt=qSeq[(m.start(0)-1):m.end(0)]
                    aaVariationList.append(AminoAcidVariation(sId,qId,pos,ref,alt,"INS",chromosome,varCount,alignInfo))
                    varCount=varCount+1;
                else: 
                    print("INSERTION at the begining of the alignment!!!");
            else:
                if qSeq[m.start(0):m.end(0)] == ("-"*span): # this is an event of deletion.
                    if m.start(0)!=0:
                        ref=sSeq[(m.start(0)-1):m.end(0)];
                        alt=qSeq[m.start(0)-1]
                        aaVariationList.append(AminoAcidVariation(sId,qId,pos,ref,alt,"DEL",chromosome,varCount,alignInfo));
                        varCount=varCount+1;
                    else: 
                        print("DELETION at the begining of the alignment!!!");
                else: ##these are asynonimous alterations or combination of INDELs and ALT.
                    ref=sSeq[m.start(0):m.end(0)];
                    alt=qSeq[m.start(0):m.end(0)];
                    ##Now these 2 strings should be checked for possible combination of INDELs and alterations.
                    if "-" not in ref: #no insertion
                        if "-" not in alt: #no deletion
                            ##all of these spaces are due to altered Amino acids
                            if len(ref)>1:
                                aaVariationList.append(AminoAcidVariation(sId,qId,pos+1,ref,alt,"ALT",chromosome,varCount,alignInfo)); # 1 is added because spaces are not representing INDELs here, so we need the start location of the event
                                varCount=varCount+1;
                            else:
                                aaVariationList.append(AminoAcidVariation(sId,qId,pos+1,ref,alt,"SAP",chromosome,varCount,alignInfo));# 1 is added because spaces are not representing INDELs here, so we need the start location of the event
                                varCount=varCount+1;
                        else: # mixture of deletion and alteration. So deletionCount does not change.                            
                            ins=re.finditer('\-+',alt);
                            prevIns=None;
                            for i in ins:
                                if i.start(0)!=0:
                                    if prevIns==None:
                                        sInd=0;
                                    else:
                                        sInd=prevIns.end(0);
                                    refPart1=ref[sInd:i.start(0)];
                                    altPart1=alt[sInd:i.start(0)];
                                    refPart2=ref[(i.start(0)-1):i.end(0)];
                                    altPart2=alt[(i.start(0)-1)];
                                    stPart2=i.start(0); #1 is not added because we want the position-1 of the current event because its event if deletion
                                    if len(refPart1)==1:
                                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+1,refPart1,altPart1,"SAP",chromosome,varCount,alignInfo));# 1 is added because these are the characters before spaces. where we need the start location of the event
                                        varCount=varCount+1;
                                    else:
                                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+1,refPart1,altPart1,"ALT",chromosome,varCount,alignInfo)); # 1 is added because these are the characters before spaces. where we need the start location of the event
                                        varCount=varCount+1;
                                    if len(refPart2)==1:
                                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+stPart2,refPart2,altPart2,"SAP",chromosome,varCount,alignInfo));
                                        varCount=varCount+1;
                                    else:
                                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+stPart2,refPart2,altPart2,"DEL",chromosome,varCount,alignInfo));
                                        varCount=varCount+1;
                                else:
                                    sp=i.end(0)-i.start(0);
                                    refPart=sSeq[m.start(0)-1:(m.start(0)+sp)]
                                    altPart=qSeq[m.start(0)-1];
                                    aaVariationList.append(AminoAcidVariation(sId,qId,pos,refPart,altPart,"DEL",chromosome,varCount,alignInfo));
                                    varCount=varCount+1;
                                prevIns=i;
                            if prevIns!=None:
                                if prevIns.end(0)<len(ref):
                                    refPart=ref[prevIns.end(0):]
                                    altPart=alt[prevIns.end(0):]
                                    aaVariationList.append(AminoAcidVariation(sId,qId,pos+1+prevIns.end(0),refPart,altPart,"ALT",chromosome,varCount,alignInfo));
                                    varCount=varCount+1;
                            else:
                                print("ERROR: This should not happen: when alt contains '-'");
                    else:
                        if "-" not in alt: #mixture of insertion and alteration. deletionCount change.
                            dels=re.finditer('\-+',ref);
                            prevDel=None;
                            intDelCount=0;
                            for i in dels:
                                if i.start(0)!=0:
                                    if prevDel==None:
                                        sInd=0;
                                    else:
                                        sInd=prevDel.end(0);
                                    refPart1=ref[sInd:i.start(0)];
                                    altPart1=alt[sInd:i.start(0)];
                                    refPart2=ref[(i.start(0)-1)];
                                    altPart2=alt[(i.start(0)-1):i.end(0)];
                                    stPart2=i.start(0)-intDelCount;
                                    if len(refPart1)==1:
                                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+1-intDelCount,refPart1,altPart1,"SAP",chromosome,varCount,alignInfo));# 1 is added because these are the characters before spaces. where we need the start location of the event
                                    else:
                                        aaVariationList.append(AminoAcidVariation(sId,qId,pos+1-intDelCount,refPart1,altPart1,"ALT",chromosome,varCount,alignInfo));# 1 is added because these are the characters before spaces. where we need the start location of the event
                                    varCount=varCount+1;
                                    aaVariationList.append(AminoAcidVariation(sId,qId,pos+stPart2,refPart2,altPart2,"INS",chromosome,varCount,alignInfo));
                                    varCount=varCount+1;
                                    deletionCount=deletionCount+len(i.group(0))
                                    intDelCount=intDelCount+len(i.group(0));
                                else:
                                    sp=i.end(0)-i.start(0);
                                    refPart=sSeq[m.start(0)-1]
                                    altPart=qSeq[m.start(0)-1:(m.start(0)+sp)];
                                    aaVariationList.append(AminoAcidVariation(sId,qId,pos,refPart,altPart,"INS",chromosome,varCount,alignInfo));
                                    varCount=varCount+1;
                                    deletionCount=deletionCount+len(i.group(0));
                                    intDelCount=intDelCount+len(i.group(0));
                                prevDel=i;
                            ##add the last one if exists, that is anything after the last "-".
                            if prevDel!=None:
                                if prevDel.end(0)<len(alt):
                                    altPart=alt[prevDel.end(0):];
                                    refPart=ref[prevDel.end(0):];
                                    aaVariationList.append(AminoAcidVariation(sId,qId,pos+1+prevDel.end(0)-intDelCount,refPart,altPart,"ALT",chromosome,varCount,alignInfo));
                                    varCount=varCount+1;
                            else:
                                print("ERROR: This should not happen: when ref contains '-'");
                        else: # mixture of alteration, insertion and deletion.
                            mixCount=mixCount+1;
                            print(qID+" mixture of insertion, deletion and alteration:"+m.group())
            prevM=m;  
        ##add the last chunk of the sequence, that comes after the last INDEL event.
        if prevM!=None: ## this should happen if there is no asynonimous ALT/SAP/INDELs, it will not be true and the whole alignment sequence should be checked for SAPs.
            qSeqPart=qSeq[prevM.end(0):];
            sSeqPart=sSeq[prevM.end(0):];
            mSeqPart=mSeq[prevM.end(0):];
            pos=prevM.end(0); ## do not subtract the deletion count as it is deducted in the synonimousSAPs function
        else:
            qSeqPart=qSeq;
            sSeqPart=sSeq;
            mSeqPart=mSeq;
            pos=0; # position is made 1 base in the synonimousSAPs function
        retRes=synonimousSAPs(qSeqPart,sSeqPart,mSeqPart,deletionCount,varCount,aaVariationList,pos,sId,qId,chromosome, alignInfo);
        aaVariationList=retRes['varList'];
        varCount=retRes['varCount'];
        #printVariationList(aaVariationList);
    else:
        print("ERROR, alignment sequences have different number of amino acids");
    return aaVariationList;

def classify(line, matchCol, evalTheshold,evalCol, gTh, gCol, lCol, qLenCol, sLenCol, qSeqCol, sSeqCol, mSeqCol, qSt, qEnd, sSt, sEnd, sId, qId, chromosome, SAPFileWriter, count, knownProteinFile, knownProteinSAPFile, isoformsFile, isoformsSAPsFile):
    mStr=line[matchCol].strip()
    #print("in classify '"+mStr+"'")
    if checkMatch(mStr,'yes')==1:
        #print("match")
        print("good match:"+str(line[gCol])+ " and threshold:"+str(gTh))
        if floatMatchCheck(float(line[gCol]),1)==1 and floatMatchCheck(float(line[lCol]),1)==1 and lengthRatioCheck(float(line[qLenCol]),float(line[sLenCol]))==1:
            print(line[qId]+" has exact match with "+line[sId]);
            knownProteinFile.write(line[qId]+","+line[sId]+"\n")
        elif evalueCheck(float(line[evalCol]),evalTheshold)==1: ##filtering out homologous alignments, generally 1*e-30
            if floatMatchCheck(float(line[gCol]),1)==1: # This means aligned seq is perfect, hence SAP within the alignment is possible, only if query/subject length is 1 or 2 AA longer.
                if floatMatchCheck(float(line[lCol]),1)==1: #this means identities/hit_length=1, here, it means that query is longer than the subject. Candidate for Isoform.
                    print(line[qId]+" is possibly an isoform of "+line[sId])
                    isoformsFile.write(line[qId]+","+line[sId]+",LongORF\n")
                elif lengthRatioCheck(float(line[qLenCol]),float(line[sLenCol]))==1:
                    print(line[qId]+" and "+line[sId]+" are of same length and partially aligned perfectly. Candidate for Isoform")
                    isoformsFile.write(line[qId]+","+line[sId]+",PartialSameLength\n")
                else:
                    ##perfect alignment but alignment is to a longer or shorter subject. Candidate for Isoform.
                    print(line[qId]+" has perfect partial alignment with "+line[sId]+", but they have different length. Candidate for Isoform or merged proteins.")
                    isoformsFile.write(line[qId]+","+line[sId]+",PartialDiffLength\n")
            elif floatAboveThreshold(float(line[gCol]),gTh)==1: ##there is SAPs or INDELs within the alignment. and if the alignment is above certain quality we check for variations.
                print("good match thresfold passed:"+line[qId])
                alignInfo={'queryStart':line[qSt],'queryEnd':line[qEnd],'qLen':line[qLenCol],'subjectStart':line[sSt],'subjectEnd':line[sEnd],'sLen':line[sLenCol],'qSerialNo':count}
                print("align Info:"+str(alignInfo))
                pId=re.search('(?:\|)(.+?)(?:\|)',line[sId]) # here group(0) will return text including '|'s, and group(1) returns only (.+?). this regex is very uniprot specific. in future, if we get reference proteome from non uniprot source, this should be revisited.
                #print("pid:"+pId.group(1))
                if pId!=None: #this should always be true here if the subject came from uniprot
                    aaVariationList=findSAPsAndINDELs(line[qSeqCol], line[sSeqCol], line[mSeqCol], pId.group(1), line[qId], line[chromosome], alignInfo)
                    printVariationListtoFile(aaVariationList,SAPFileWriter)
                else:
                    print(line[qId]+" mapped to "+line[sId]+" which did not have uniprot style id.")
                    aaVariationList=findSAPsAndINDELs(line[qSeqCol], line[sSeqCol], line[mSeqCol], line[sId], line[qId], line[chromosome], alignInfo)
                    printVariationListtoFile(aaVariationList,SAPFileWriter)
                if int(line[qSt])==1 and int(line[sSt])==1 and int(line[qEnd])==int(line[qLenCol]) and int(line[sEnd])==int(line[sLenCol]): #equality of the query and the subject is not enough to say that the alignment covers the whole sequences hence checking that end points of the alignments are the positions of the sequences.
                    print("knownProteinSap")
                    knownProteinSAPFile.write(line[qId]+","+line[sId]+"\n")
                else:
                    if floatMatchCheck(float(line[lCol]),1)==1:
                        isoformsSAPsFile.write(line[qId]+","+line[sId]+",LongORF\n");
                    elif lengthRatioCheck(float(line[qLenCol]),float(line[sLenCol]))==1:
                        isoformsSAPsFile.write(line[qId]+","+line[sId]+",SameLength\n");
                    else:
                        isoformsSAPsFile.write(line[qId]+","+line[sId]+",DiffLength\n");
            else:
                print(line[qId]+" has poor alignment quality, quality:"+line[gCol]);
        else:
            print(line[qId]+" did not pass e-value threshold.")

def read(filename, matchCol, evalTheshold, evalCol, gTh, gCol, lCol, qLenCol, sLenCol, qSeqCol, sSeqCol, mSeqCol, qSt, qEnd, sSt, sEnd, sId, qId, chromosome, SAPFileName, knownProteinFileName, knownProteinSAPFileName, isoformsFileName, isoformsSAPsFileName):
    with open(filename, newline='') as csvfile, open(SAPFileName, 'w', newline='') as SAPFile, open(knownProteinFileName, 'w', newline='') as knownProteinFile, open(knownProteinSAPFileName, 'w', newline='') as knownProteinSAPFile, open(isoformsFileName, 'w', newline='') as isoformsFile, open(isoformsSAPsFileName, 'w', newline='') as isoformsSAPsFile:
        reader = csv.reader(csvfile, delimiter=',')#,quoting=csv.QUOTE_NONNUMERIC
        #SAPFileWriter = open(SAPFile, delimiter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        count=0
        columns=[]
        SAPFile.write(AminoAcidVariation.printHeader()+"\n");
        knownProteinFile.write("ORF Id, Protein ID\n");
        knownProteinSAPFile.write("ORF Id, Protein ID\n");
        isoformsFile.write("ORF Id, Protein ID, Type\n");
        isoformsSAPsFile.write("ORF Id, Protein ID, Type\n");
        for line in reader:
            if count>0:
                #print("classify called");
                classify(line, matchCol, evalTheshold,evalCol, gTh, gCol, lCol, qLenCol, sLenCol, qSeqCol, sSeqCol, mSeqCol, qSt, qEnd, sSt, sEnd, sId, qId, chromosome, SAPFile, count, knownProteinFile, knownProteinSAPFile, isoformsFile, isoformsSAPsFile)
                count=count+1
            else:
                count=count+1

filename="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_with_LocationV2.csv"
SAPFileName="human_adeno_mydb_pasa.assemblies_ORFs_with_Location_VariationV5.vcf"
knownProteinFileName="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsV5.csv"
knownProteinSAPFileName="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsSAPsV5.csv"
isoformsFileName="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsV5.csv"
isoformsSAPsFileName="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsSAPsV5.csv"
matchCol=2
evalTheshold=0.000000000000000000000000000001
evalCol=6
gTh=0.5
gCol=11
lCol=12
qLenCol=1
sLenCol=5
qSeqCol=19
sSeqCol=20
mSeqCol=21
qSt=15
qEnd=16
sSt=17
sEnd=18
sId=4
qId=0
chromosome=22
alignCol=9
read(filename, matchCol, evalTheshold, evalCol, gTh, gCol, lCol, qLenCol, sLenCol, qSeqCol, sSeqCol, mSeqCol, qSt, qEnd, sSt, sEnd, sId, qId, chromosome, SAPFileName, knownProteinFileName, knownProteinSAPFileName, isoformsFileName, isoformsSAPsFileName)


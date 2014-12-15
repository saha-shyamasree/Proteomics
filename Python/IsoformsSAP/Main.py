import os
import csv
import re

def shortORFFileWrite(line, shortORFFileName):
    with open(shortORFFileName, 'a', newline='') as shortORFFile:
        shortORFFileWriter = csv.writer(shortORFFile, delimeter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        shortORFFileWriter.writerow(line)
        shortORFFile.close();

def extendedORFFileWrite(line, extendedFileName):
    with open(extendedORFFileName, 'a', newline='') as extendedORFFile:
        extendedORFFileWriter = csv.writer(extendedORFFile, delimeter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        extendedORFFileWriter.writerow(line)
        extendedORFFile.close()

def variationLengthORFWrite(line, qSeqCol, sSeqCol, extendedORFFileName, shortORFFileName):
    if len(line[qSeqCol])>len(line[sSeqCol]):
        extendedORFFileWrite(line, extendedFileName)
    else:
        shortORFFileWrite(line, shortORFFileName)

def SAPORFWrite(line, SAPFileName):
    with open(SAPFileName, 'a', newline='') as SAPFile:
        SAPORFFileWriter = csv.writer(SAPFile, delimeter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        SAPORFFileWriter.writerow(line)
        SAPORFFile.close()

def unidentified(line, unidentifiedFileName):
    with open(unidentifiedFileName, 'a', newline='') as unidentifiedFile:
        unidentifiedFileWriter = csv.writer(unidentifiedFile, delimeter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        unidentifiedFileWriter.writerow(line)
        unidentifiedFile.close()

def lengthBaseSAPCheck(line, qSeqCol, sSeqCol, qSeqLengthCol, qSeqLengthCol, extendedORFFileName, shortORFFileName, unidentifiedFileName):
    sSeq=line[sSeqCol] #[(sStart-1):sEnd]  ##because first index inclusive and last index exclusive and python index is 0 based wheras blast index is 1 based
    qSeq=line[qSeqCol] #[(qStart-1):qEnd]
    
    if (line[qSeqLengthCol]==line[qSeqLengthCol]+1 and line[qSeqLengthCol]==len(qSeq)) or (line[sSeqLengthCol]+1==line[qSeqLengthCol] and line[sSeqLengthCol]==len(sSeq)):
        SAPORFWrite(line, SAPFileName)
    elif line[sSeqLengthCol]==line[qSeqLengthCol]+2:
        if line[qSeqLengthCol]==len(qSeq):
            if sStart==2:  #if length of the subject sequence is 2 amino acids longer than the query sequence, and the alignment is equal to query sequence, plus alignment starts from second amino acid that means 2 extra amino acids are in 2 ends.
                SAPORFWrite(line, SAPFileName)
            else:
                shortORFFileWrite(line, shortORFFileName)
        else:
            unidentified(line, unidentifiedFileName)
    elif line[qSeqLengthCol]==line[sSeqLengthCol]+2:
        if line[sSeqLengthCol]==len(sSeq):
            if qStart==2:
                SAPORFWrite(line, SAPFileName)
            else:
                extendedORFFileName(line, extendedORFFileName)
        else:
            unidentified(line, unidentifiedFileName)
    elif line[qSeqLengthCol]==line[sSeqLengthCol]:  ###that is query and subject is of same length, but smaller part of it match.
        if line[qSeqLengthCol]==len(qSeq)+1 or line[sSeqLengthCol]==len(sSeq)+1:
            SAPORFWrite(line, SAPFileName)
    elif line[qSeqLengthCol]==len(qSeq) or line[sSeqLengthCol]==len(sSeq):
        variationLengthORFWrite(line, qSeqCol, sSeqCol, extendedORFFileName, shortORFFileName)

def characterCheckLoop(seq1, seq2):
    if len(seq1)==len(seq2): ##should be equal size
        SapGapMap = {}
        st=-1
        end=-1
        for i in range(0,len(seq1):
            if seq1[i]!=seq2[i]:
                if st==-1:
                    st=i
                    end=i
                else:
                    end=i
            else:
                #if st is not -1 that means there was mismatch before this, so that should go into SapGapMap.
                st=-1 #???
    else:
        print("ERROR")

def checkCharacterByCharacterCheck(line, qSeqCol, sSeqCol, gapCol, SAPFileName):
    sSeq=line[sSeqCol]#[(sStart-1):sEnd]  ##because first index inclusive and last index exclusive and python index is 0 based wheras blast index is 1 based
    qSeq=line[qSeqCol]#[(qStart-1):qEnd]
    #use re to find locations of gaps, for each gap count the number of '-' character for each group.
    pattern=re.compile("-+")
    prevEnd=-1
    iterator=pattern.finditer(line[qSeqCol])
    if line[gapCol]==1:
        indx=line[qSeqCol].find("-")
        if indx==-1:
            indx=line[sSeqCol].find("-")
            if indx!=-1:
                if line[qSeqCol][0:indx]==line[sSeqCol][0:indx]:
                    if line[qSeqCol][indx+1:]==line[sSeqCol][indx+1:]:
                        SAPORFWrite(line, SAPFileName)
                    else:
                        print("character by character check")
                        #characterCheckLoop(line[qSeqCol][indx+1:], line[sSeqCol][indx+1:])
                else:
                    print("character by character check")
                    #characterCheckLoop(line[qSeqCol][0:indx], line[sSeqCol][0:indx])
            else:
                print("ERROR, wrong gap count")
        else:
            if line[qSeqCol][0:indx]==line[sSeqCol][0:indx]:
                if line[qSeqCol][indx+1:]==line[sSeqCol][indx+1:]:
                    SAPORFWrite(line, SAPFileName)
                else:
                    print("character by character check")
                    characterCheckLoop(line[qSeqCol][indx+1:], line[sSeqCol][indx+1:])
            else:
                print("character by character check")
                characterCheckLoop(line[qSeqCol][0:indx], line[sSeqCol][0:indx])
                if line[qSeqCol][indx+1:]!=line[sSeqCol][indx+1:]:
                    characterCheckLoop(line[qSeqCol][indx+1:], line[sSeqCol][indx+1:])
                    ##based on the return value from characterCheckLoop, apropriate function to write the line should be called.
    elif line[gapCol]>1:
        
        for p in pattern.finditer(line[qSeqCol]):
            tempSpan=p.span() ##collect these locations. these spans give us length of each gap and allow us to compare segments divided by these gaps. Same comparisons shoud take place for subject sequence. we can check whether subject quesry has any gap by counting number of characters in those patterns.
            if prevEnd==-1:
                #print("character by character check")
                if line[qSeqCol][0:tempSpan[0]]!=line[qSeqCol][0:tempSpan[0]]:
                    characterCheckLoop(line[qSeqCol][0:tempSpan[0]],line[sSeqCol][0:tempSpan[0]])
            else:
                print("character by character check")
                if line[qSeqCol][prevEnd:tempSpan[0]]!=line[qSeqCol][prevEnd:tempSpan[0]]:
                    characterCheckLoop(line[qSeqCol][0:tempSpan[0]],line[sSeqCol][0:tempSpan[0]])
            prevEnd=tempSpan[1]
    elif line[gapCol]==0:
        characterCheckLoop(line[qSeqCol], line[sSeqCol])
    
def classify(line, evalTheshold,evalCol, gTh, gCol, lCol, qLenCol, sLenCol, gapCol, qSeqCol, sSeqCol, exactMatchFilename, extendedORFFileName, shortORFFileName):
    if evalueCheck(line[evalCol],evalTheshold)==1:
        if goodMatchCheck(line[gCol],1)==1 and longMatchCheck(line[lCol])==1 and lengthRatioCheck(qLenCol,sLenCol)==1:
            exactMatchFileWriter(line, exactMatchFilename)
        elif goodMatchCheck(line[gCol],gTh)==1:
            sSeq=line[sSeqCol] #[(sStart-1):sEnd]  ##because first index inclusive and last index exclusive and python index is 0 based wheras blast index is 1 based
            qSeq=line[qSeqCol] #[(qStart-1):qEnd]
            if sSeq == qseq:
                lengthBaseSAPCheck(line, qSeqCol, sSeqCol, extendedORFFileName, shortORFFileName)
            elif longMatchCheck(line[lCol],1) == 1:
                
            else:
                checkCharacterByCharacterCheck(line, qSeqCol, sSeqCol, gapCol)

def read(filename,evalTheshold,evalCol, gTh, gCol, lCol, qLenCol, sLenCol, qSeqCol, sSeqCol, exactMatchFilename, extendedORFFileName, shortORFFileName):
    with open(filename, newline='') as csvfile, open(exactMatchFilename, 'w', newline='') as exactMatchFile, open(extendedORFFileName, 'w', newline='') as extendedORFFile, open(shortORFFileName, 'w', newline='') as shortORFFile:
        reader = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        exactMatchFileWriter = csv.writer(exactMatchFile, delimeter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        extendedORFFileWriter = csv.writer(extendedORFFile, delimeter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        shortORFFileWriter = csv.writer(shortORFFile, delimeter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        exactMatchFile.close()
        extendedORFFile.close()
        shortORFFile.close()
        
        count=0
        for line in reader:
            if count>0:
                classify(line, evalTheshold,evalCol, gTh, gCol, lCol, qLenCol, sLenCol, exactMatchFilename, extendedORFFileName, shortORFFileName)
            else:
                count=count+1
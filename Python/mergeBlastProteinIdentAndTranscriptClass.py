## this code takes list of ORFs grouped based of the biotype of their parent transcripts, ORFs ids mapped to uniprot using blast, protein identification results, and SAP/INDELs to merge these information. #Output of this script will then be merged with peptide evidences.
import csv
import re

## column numbers of necessary fields in different files.
prtGrpCol=2
orfSubjectCol=
orfLengthRatioCol=
transcriptSubjectCol=
transcriptLengthRatioCol=
locationCol=


def compareProteinParentTranscriptClass(prtSubGrpFilepath, orfToUniprotMap, prtIdCol):
    with open(filePath, newline='') as transClassFile:
        count=0
        
        for line in transClassReader:
            if count>0:
                ##find the ORF in orfToUniprotMap
                
        

def mergeInfo(ORFsBioTypeProtein, ORFsBioTypePeptide, ORFstoProteomeBlast, SAPs):
    ORFsBioTypeProteinCSV=csv.reader(ORFsBioTypeProtein, delimiter="\t")
    ORFsBioTypePeptideCSV=csv.reader(ORFsBioTypePeptide, delimiter="\t")
    ORFstoProteomeBlastCSV=csv.reader(ORFstoProteomeBlast, delimiter="\t" )
    SAPsCSV=csv.reader(SAPs, delimiter="\t")
    for line in ORFsBioTypeProteinCSV:
        ##find uniprot mapping and relavent informations, peptides, and any SAPs
        
    
    

def readFiles(ORFsBioTypeProteinFilePath, ORFsBioTypePeptideFilePath, ORFstoProteomeBlastFilePath, SAPFilePath):
    with open(ORFsBioTypeProteinFilePath, 'r') as ORFsBioTypeProteinFile, open(ORFsBioTypePeptideFilePath, 'r') as ORFsBioTypePeptideFile, open(ORFstoProteomeBlastFilePath, 'r') as ORFstoProteomeBlastFile, open(SAPFilePath, 'r') as SAPFile:
        ORFsBioTypeProtein = ORFsBioTypeProteinFile.readlines()
        ORFsBioTypePeptide = ORFsBioTypePeptideFile.readlines()
        ORFstoProteomeBlast = ORFstoProteomeBlastFile.readlines()
        SAPs = SAPFile.readlines()
        mergeInfo(ORFsBioTypeProtein, ORFsBioTypePeptide, ORFstoProteomeBlast, SAPs)

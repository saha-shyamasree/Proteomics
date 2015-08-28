## This Python code Reads (1) the list of known proteins with/without SAPs, (2) isoforms with/without SAPs, and (3) the list of ORFs classified according to their
## transcript bio-type. Aim of this code is to exclude (1) from the list of (3).

import pandas as pd
import re

def readList(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep)
    return fileDFObj;

def replaceCharacter(Mat, column, character, replace):
    ## this function replace characters like ',',';' with replace character.
    Mat[column]=Mat[column].str.replace(character, replace)
    return Mat[column]

def subsetProteins(Mat1, MatList, column1, column2):
    ## this code find ORFs from MatList, that are in Mat1 and returns MatList[i]-(MatList[i] Intersection Mat1)
    subMats=[]
    subCount=[]
    for i in range(len(MatList)):
        ## count tells us how many of MatList[i] are in Mat1
        subCount[i]=MatList[i][column2].isin(Mat1[column1]).sum()
        SubMats[i]=MatList[i][~MatList[i][column2].isin(Mat1[column1])]
    
    return {'Count':subCount,'Matrix':subMats}
    
    
def main(knownPro, knownProSAP, iso, isoSAP, nonCDSFileNameList):
    ## Read a) known proteins, b) known proteins with SAPs including ALT/INDELs, c) isoforms, and d) isoforms with SAPs/ALT/INDELs
    ## Known Protein
    knownProteins=replaceCharacter(readList(knownPro,','),';',',')
    ## Known with SAPs/ALT/INDELs
    knownProteinsSAPs=replaceCharacter(readList(knownProSAP,','),';',',')
    ## Isoforms
    isoforms=replaceCharacter(readList(iso,','),';',',')
    ## Isoforms with SAPs/ALT/INDELs
    isoformsSAPs=replaceCharacter(readList(isoSAP,','),';',',')
    
    ## Read i) anisense, ii) linc, iii) nonsense-mediated decay, iv) processed transcript, v) reatined intron and any other ORFs from
    ## non CDS RNA bio type.
    nonCDS=[]
    for i in range(length(nonCDSFilenameList)):
        nonCDS[i]=readList(nonCDSFilenameList[i],'\t')
    
    ##MSGF+ protein file has comma in protein accession, whereas MatList files, i.e. known proteins etc. have ';'.
    
    ## Remove all proteins/ORFs from i-v if they are found in a. Make separate lists of those ORFs/proteins from the i-v lists that
    ## exist in b, c or d.
    ##
    res1=subsetProteins(knownProteins,nonCDS)
    res2=subsetProteins(knownProteinsSAPs,res1['Matrix'])
    res3=subsetProteins(isoforms,res2['Matrix'])
    res4=subsetProteins(isoformsSAPs,res3['Matrix'])

knownPro="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsV7.csv"
knownProSAP="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsSAPsV7.csv"
iso="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsV7.csv"
isoSAP="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsSAPsV7.csv"

nonCDSPrtFileNameList=["D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonsense_mediated_decayPrt.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/processed_transcriptPrt.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/retained_intronPrt.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/antisensePrt.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/lincPrt.tsv"]
nonCDSPrtFileNameList=["D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonsense_mediated_decay.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/processed_transcript.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/retained_intron.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/antisense.tsv","D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/linc.tsv"]

###This code can count number of unique peptides for a list of protein groups. Input is in csv from mzidentml export type protein
import csv

def readProteinGroup(filePath):
    with open(filePath, 'r', newline='') as RD:
        csv_rd=csv.DictReader(RD, delimiter="\t")
        peptides_list=[]
        for line in csv_rd:
            if line['unique.peptides']:
                peptides_list=peptides_list+line['unique.peptides'].split(';')
    return peptides_list;

##CDS
filePath="D:/data/Results/Human-Adeno/Identification/PASA/sORF/CDSPrt.tsv"
peptides=readProteinGroup(filePath)
print("CDS:"+str(len(peptides)))

##nonsense mediated decay
filePath="D:/data/Results/Human-Adeno/Identification/PASA/sORF/nonsense_mediated_decayPrt.tsv"
peptides=readProteinGroup(filePath)
print("nonsense mediated decay:"+str(len(peptides)))

##processed transcripts
filePath="D:/data/Results/Human-Adeno/Identification/PASA/sORF/processed_transcriptPrt.tsv"
peptides=readProteinGroup(filePath)
print("processed transcript:"+str(len(peptides)))

##retained intron
filePath="D:/data/Results/Human-Adeno/Identification/PASA/sORF/retained_intronPrt.tsv"
peptides=readProteinGroup(filePath)
print("retained intron:"+str(len(peptides)))

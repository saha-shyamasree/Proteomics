### GFF3 comparison
import sys
import pandas as pd
import argparse
import numpy as np

def readFile(filename, sep, headerFlag):
    if headerFlag==0:
        fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    elif headerFlag==1:
        fileDFObj = pd.read_table(filename, sep=sep, header=None, keep_default_na=False, na_values=[''])
    else:
        print("Unrecognized Header Flag")
    return fileDFObj;

#AFile=Sample gene structure, gff3 format
#BFile=Current gene annotation of the species, gff3 format
def compareGFF(AFile,BFile):
    ##Create Bprime, i.e. take the standard annotatioon file and remove all the features that has no overlap with the sample gene structure.
    command1="intersectBed -u -s -a "+BFile+" -b "+AFile+" > "+BFile+"_Overlap.gff3"
    ##A subtract B gives us that A has but not B.
    command2="subtractBed -s -a "+AFile+" -b "+BFile+" > "+AFile+"_Extended.gff3"
    ##B' subtract A gives us that B' has but not A.
    command3="subtractBed -s -a "+BFile+"_Overlap.gff3"+" -b "+AFile+" > "+AFile+"_Missing.gff3"
    
    a-b=AFile+"_Extended.gff3"
    bprime-a=AFile+"_Missing.gff3"

def featureAlignment(pdObjGrp, idCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol):
    featureAClass=['gene','mRNA','CDS','exon','three_prime_utr','five_prime_utr']
    featureBClass=['gene','mRNA','CDS','exon','UTR','start_codon']
for name, alignments in pdObjGrp:
    genes=alignments[alignments.iloc[:,aFeatureCol]=='gene']
    if len(alignments)>1:
        #print(alignments.iloc[:,aStartCol])
        starts=np.array([max(alignments.iloc[i,aStartCol],alignments.iloc[i,bStartCol]) for i in range(0,len(alignments))])
        stops=np.array([min(alignments.iloc[i,aStopCol],alignments.iloc[i,bStopCol]) for i in range(0,len(alignments))])
        overlaps=stops-starts
        lenRatioDis=1-((0.5*(overlaps/(alignments.iloc[:,aStopCol]-alignments.iloc[:,aStartCol])))+(0.5*(overlaps/(alignments.iloc[:,bStopCol]-alignments.iloc[:,bStartCol]))))
        #for a in alignments:
        exact=alignments[lenRatioDis==0]
        if exact.empty:
            print(name+": no exact overlap")
        else:
            print(exact.size)
        ### Filtering alignments that has same or similar feature type.
        #print(alignments)

for name, alignments in genesGrp:
    #genes=alignments[alignments.iloc[:,aFeatureCol]=='gene']
    if len(alignments)>1:
        #print(alignments)
        #print(alignments.iloc[:,aStartCol])
        starts=np.array([max(alignments.iloc[i,aStartCol],alignments.iloc[i,bStartCol]) for i in range(0,len(alignments))])
        stops=np.array([min(alignments.iloc[i,aStopCol],alignments.iloc[i,bStopCol]) for i in range(0,len(alignments))])
        overlaps=stops-starts
        lenRatioDis=1-((0.5*(overlaps/(alignments.iloc[:,aStopCol]-alignments.iloc[:,aStartCol])))+(0.5*(overlaps/(alignments.iloc[:,bStopCol]-alignments.iloc[:,bStartCol]))))
        #for a in alignments:
        print(lenRatioDis)
        exact=alignments[lenRatioDis==0]
        print(exact)
        ### Filtering alignments that has same or similar feature type.
        #alignments

filename="D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.genome.gff3_identified_intersect_trans_ref_wa_wb.gff3"
ann=readFile(filename, '\t', 1)
idCol=8
aStartCol=3
aStopCol=4
aFeatureCol=2
bStartCol=12
bStopCol=13
bFeatureCol=11
pdObjGrp=ann.groupby(8)
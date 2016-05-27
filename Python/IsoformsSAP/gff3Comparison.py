### GFF3 comparison
import sys
import pandas as pd
import argparse
import numpy as np
import re

class Gene:
    """This clas holds a gene structure. The purpose of this class is to hold all the overlaps of a gene from dataset B."""
    def __init__(self,gid, alignments, mRNAs, exons, cds, utr_three, utr_five):
        self.gid=gid
        self.alignments=alignments
        self.mRNAs=mRNAs
        self.exons=exons
        self.cds=cds
        self.utr_three=utr_three
        self.utr_five=utr_five
    def getGeneID(self):
        return self.gid
    def getAlignments(self):
        return self.alignments
    def getmRNAs(self):
        return self.mRNAs
    def getExons(self):
        return self.exons
    def getCDS(self):
        return self.cds
    def getUTRThree(self):
        return self.utr_three
    def getUTRFive(self):
        return self.utr_five
    def setGeneID(self,gid):
        self.gid=gid
    def setAlignments(self,alignments):
        self.alignments=alignments
    def setmRNAs(self,mRNAs):
        self.mRNAs=mRNAs
    def setExons(self,exons):
        self.exons=exons
    def setCDS(self, cds):
        self.cds=cds
    def setUTRThree(self, utr_three):
        self.utr_three=utr_three
    def setUTRFive(self, utr_five):
        self.utr_five=utr_five

    
class mRNA:
    """This clas holds a mRNA structure. The purpose of this class is to hold all the overlaps of a mRNA from dataset B."""
    def __init__(self,mid, alignments, parent, exons, cds, utr_three, utr_five):
        self.mid=mid
        self.alignments=alignments
        self.parent=parent
        self.exons=exons
        self.cds=cds
        self.utr_three=utr_three
        self.utr_five=utr_five
    def getmRNAID(self):
        return self.mid
    def getAlignments(self):
        return self.alignments
    def getParent(self):
        return self.parent
    def getExons(self):
        return self.exons
    def getCDS(self):
        return self.cds
    def getUTRThree(self):
        return self.utr_three
    def getUTRFive(self):
        return self.utr_five
    def setGeneID(self,gid):
        self.gid=gid
    def setAlignments(self,alignments):
        self.alignments=alignments
    def setParent(self,parent):
        self.parent=parent
    def setExons(self,exons):
        self.exons=exons
    def setCDS(self, cds):
        self.cds=cds
    def setUTRThree(self, utr_three):
        self.utr_three=utr_three
    def setUTRFive(self, utr_five):
        self.utr_five=utr_five

class SubFeature:
    """This clas holds sub-features such as exons, cds, UTRs. The purpose of this class is to hold all the alignments of a sub-feature from dataset B."""
    def __init__(self, fid, alignments, parent, featureType):
        self.fid=fid
        self.alignments=alignments
        self.parent=parent
        self.featureType=featureType
    def getFeatureID(self):
        return self.fid
    def getAlignments(self):
        return self.alignments
    def getParent(self):
        return self.parent
    def getFeatureType(self):
        return self.featureType

    

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
'''
def compareGFF(AFile,BFile):
    ##Create Bprime, i.e. take the standard annotatioon file and remove all the features that has no overlap with the sample gene structure.
    command1="intersectBed -u -s -a "+BFile+" -b "+AFile+" > "+BFile+"_Overlap.gff3"
    ##A subtract B gives us that A has but not B.
    command2="subtractBed -s -a "+AFile+" -b "+BFile+" > "+AFile+"_Extended.gff3"
    ##B' subtract A gives us that B' has but not A.
    command3="subtractBed -s -a "+BFile+"_Overlap.gff3"+" -b "+AFile+" > "+AFile+"_Missing.gff3"
    
    a-b=AFile+"_Extended.gff3"
    bprime-a=AFile+"_Missing.gff3"
'''

def checkChildAnnotation(name, BMaps, aIDCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap, maps):
    ##if name is a gene, we can have direct relation with mRNAs only. For each mRNA, we need to find the children.
    match=re.search("(?:ID=)([^;]*)",name)
    if match and len(match.groups())>0:
        featureID=match.group(1)
        ##look for children
        if len(BMaps)>0:
            if BMaps[aFeatureCol].iloc[0]=='gene':
                ##this is a gene
                sStr=re.escape("Parent="+featureID+";")
                children=maps[maps[aIDCol].str.contains(sStr)]
                cGrp=children.groupby(aIDCol)
                #finding alignemnt of only these mRNAs might not be enough.
                selectMap(cGrp, aIDCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap, maps)
            elif BMaps[aFeatureCol].iloc[0]=='mRNA':
                ##this is a mRNA
                sStr=re.escape("Parent="+featureID)+"$"
                children=maps[maps[aIDCol].str.contains(sStr)]
                cGrp=children.groupby(aIDCol)
                #finding alignemnt of only these mRNAs might not be enough.
                selectMap(cGrp, aIDCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap, maps)
            #


def smalletLengthRatio(unqLenRatio, delta):
    selectedLengthRatio=[unqLenRatio[0]]
    for i in range(1,len(unqLenRatio)):
        if unqLenRatio[i]-unqLenRatio[i-1]<=delta:
            ##Difference between two consecutive score is less than delta. hence we need to keep this length ratio.
            selectedLengthRatio.append(unqLenRatio[i])
        else:
            break
    return selectedLengthRatio
    

def findMapsToConsider(alignments, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta):
    starts=np.array([max(alignments.iloc[i,aStartCol],alignments.iloc[i,bStartCol]) for i in range(0,len(alignments))])
    stops=np.array([min(alignments.iloc[i,aStopCol],alignments.iloc[i,bStopCol]) for i in range(0,len(alignments))])
    overlaps=stops-starts
    lenRatioDis=1-((0.5*(overlaps/(alignments.iloc[:,aStopCol]-alignments.iloc[:,aStartCol])))+(0.5*(overlaps/(alignments.iloc[:,bStopCol]-alignments.iloc[:,bStartCol]))))
    lenRatioDis.sort()
    unqLengthRatioDis = lenRatioDis.unique().tolist()
    ##Only considering the minRatio might not be enough, speacilly for long features. I need to have delta that checks the difference bettween minRation and the secondMinRatio.
    ##It the difference is less than the delta, its not good enough to consider only the minRatio ones.
    ##Delta is static now, but in future it should be calculated based on the length ratio difference distribution.
    selectedRatio=smalletLengthRatio(unqLengthRatioDis, delta)
    ##we want the max from the selected length ratios.
    maxLenRatio=max(selectedRatio)
    print("max LenRatio:"+str(maxLenRatio))
    #print(lenRatioDis[lenRatioDis==minLenRatio])
    indices=lenRatioDis[lenRatioDis<=maxLenRatio].index.tolist()
    return indices
    
def selectMap(grpObj,idCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap, parentMaps):
    count1=0
    countMore=0
    for name, alignments in grpObj:
        #genes=alignments[alignments.iloc[:,aFeatureCol]=='gene']
        if len(alignments)>1:
            #print(alignments.iloc[:,aStartCol])
            
            exact=findMapsToConsider(alignments, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta)
            if len(exact)<=0:
                print(name+": Error 1, at least one map should be there")
            else:
                print("Total alignment:"+str(len(alignments)))
                print("Number of low distant maps: "+str(len(exact)))
                if len(exact)==1:
                    ##this solves the problem of chosing.
                    print("One winner")
                    countMore=countMore+1
                    featureMap=featureMap.append(alignments.loc[exact],ignore_index=True)
                    print("Length of featureMap:"+str(len(featureMap)))
                else:
                    if alignments[aFeatureCol].iloc[0]=='gene':
                        print("Gene")
                        ##this means this is a gene. Hence, our decision will ba based on the feature type of B.
                        subAlignments=alignments.loc[exact]
                        print("Total sub alignment:"+str(len(subAlignments)))
                        geneMaps=subAlignments[subAlignments[bFeatureCol].str.contains('gene')]
                        print("Number of gene maps:"+str(len(geneMaps)))
                        if len(geneMaps)==0:
                            ##This means none of the B features is gene.
                            ##Now check for next feature type, which is mRNA
                            rnaMaps=subAlignments[subAlignments[bFeatureCol].str.contains('mRNA')]
                            if len(rnaMaps)==0:
                                print("PROBLEM, gene should atleast map to one mRNA")
                            elif len(rnaMaps)==1:
                                featureMap=featureMap.append(rnaMaps,ignore_index=True)
                            else:
                                ##check how children features are mapped to the features from B. 
                                checkChildAnnotation(name, rnaMaps, aIDCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap, maps)
                        elif len(geneMaps)==1:
                            featureMap=featureMap.append(geneMaps,ignore_index=True)
                        else:
                            print(name +" : Thinking TIME")
                            print(subAlignments)
                            checkChildAnnotation(name, geneMaps, aIDCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap, maps)
                    elif alignments[aFeatureCol].iloc[0]=='mRNA':
                        subAlignments=alignments.loc[exact]
                        mRNAMaps=subAlignments[subAlignments[bFeatureCol].str.contains('mRNA')]
                        print("Number of mRNA maps:"+str(len(mRNAMaps)))
                        if len(mRNAMaps)==0:
                            ##this mRNA did not map to any mRNA, check whether it maps to gene
                            print("Count how many times this is happening")
                            '''
                            geneMaps=subAlignments[subAlignments[bFeatureCol].str.contains('gene')]
                            if len(geneMaps)==0:
                                print("PROBLEM, mRNA dd not map to any mRNA, "+name+" does not have a map to a gene either")
                            elif len(geneMaps)==1:
                                print("gene map for the mRNA was selected")
                                featureMap=featureMap.append(geneMaps,ignore_index=True)
                            else:
                               
                            '''
                        elif len(mRNAMaps)==1:
                            featureMap=featureMap.append(mRNAMaps,ignore_index=True)
                        else:
                            checkChildAnnotation(name, rnaMaps, aIDCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap, maps)
                    elif alignments[aFeatureCol].iloc[0]=='exon':
                        ##this does not have a child
                        subAlignments=alignments.loc[exact]
                        print("Total sub alignment (exon):"+str(len(subAlignments)))
                        ##keep only those alignments which has same parent B feature as its parent. i.e. if parent of this exon is mRNA1 and mRNA1 has mapped to ENS1, ENS2
                        ##then for this exon if alignemnts are exon1.ENS1, exon1.ENS2, exon1.ENS3, keep only exon1.ENS1, exon1.ENS2.
                        parentsFeatureMap=re.escape(parentMaps[bIDCol].str.extract("(?:ID=)([^;]*)")) ### need to extract the ID=.... part. Check whether its working
                        parentStr="|".join(parentsFeatureMap)
                        subAlignmentsParent=subAlignments[subAlignments[bIDCol].str.contains(parentStr)]
                        ###I should return this, so that I can count number of sub-feature per gene/mRNA and decide the final map.
                        ##Some  exon/CDS/three_prime/five_prime might not have any parent in the interset file. I need to check those.
                    elif alignments[aFeatureCol].iloc[0]=='CDS':##am I treating exon/cds/3'/5' diffecrently?? if not I should not create separate if.
                ### Check whether the parent feature has been added to the featureMap. 
        else:
            if len(alignments)==1:
                count1=count1+1
                featureMap=featureMap.append(alignments.iloc[0],ignore_index=True)
            else:
                print(name+": Error 2, at least one map should be there")
    print("Count1:"+str(count1))
    print("CountMore:"+str(countMore))
    return featureMap



'''
def featureAlignment(ann, pdObjGrp, idCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, newAnnotation):
    featureAClass=['gene','mRNA','CDS','exon','three_prime_utr','five_prime_utr']
    featureBClass=['gene','mRNA','CDS','exon','UTR','start_codon']
    ## Only gene or mRNA can be parent.

print("number of features:"+str(len(pdObjGrp)))

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
'''
filename="D:\data\Oliver\PASA\G10\G10.assemblies.fasta.transdecoder.genome.gff3_identified_intersect_trans_ref_wa_wb.gff3"
maps=readFile(filename, '\t', 1)
aIDCol=8
aStartCol=3
aStopCol=4
aFeatureCol=2
bStartCol=12
bStopCol=13
bFeatureCol=11
bIDCol=17
delta=0.01
pdObjGrp=maps.groupby(8)
genes=maps[maps[aFeatureCol].str.contains('gene')]
mRNA=maps[maps[aFeatureCol].str.contains('mRNA')]
    
count1=0
countMore=0
featureMap=pd.DataFrame();
gGrp=genes.groupby(aIDCol)
mGrp=mRNA.groupby(aIDCol)

selectMap(gGrp,aIDCol, aStartCol, aStopCol, aFeatureCol, bStartCol, bStopCol, bFeatureCol, delta, featureMap)
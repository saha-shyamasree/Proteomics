## This Python code Reads (1) the list of known proteins with/without SAPs, (2) isoforms with/without SAPs, and (3) the list of ORFs classified according to their
## transcript bio-type. Aim of this code is to exclude (1) from the list of (3).

import pandas as pd
import re
import os

def readList(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep)
    return fileDFObj;

def replaceCharacter(Mat, column, character, replace):
    ## this function replace characters like ',',';' with replace character.
    Mat[column]=Mat[column].str.replace(character, replace)
    #print(Mat[column].head(5))
    return Mat

def subsetProteins(Mat1, MatList, column1):
    ## this code find ORFs from MatList, that are in Mat1 and returns MatList[i]-(MatList[i] Intersection Mat1)
    subMats=[None]*len(MatList)
    subCount=[None]*len(MatList)
    for i in range(len(MatList)):
        ## count tells us how many of MatList[i] are in Mat1
        #print(str(i))
        Mat2=MatList[i]
        #colnames=Mat2.columns.values
        #print(colnames)
        subCount[i]=Mat2['description'].isin(Mat1[column1]).sum()
        subMats[i]=Mat2[~Mat2['description'].isin(Mat1[column1])]
        print("Main ORF:"+str(subMats[i]['description'].str.contains("Dataset_A").sum()))
        print("3'uORF:"+str(subMats[i]['description'].str.contains("Dataset_C").sum()))
        print("5'uORF:"+str(subMats[i]['description'].str.contains("Dataset_B").sum()))
    return {'Count':subCount,'Matrix':subMats}

def subsetPSMs(Mat1List, Mat2List, column1):
    ## this code find ORFs from MatList, that are in Mat1 and returns MatList[i]-(MatList[i] Intersection Mat1)
    subMats=[None]*len(Mat1List)
    subCount=[None]*len(Mat1List)
    for i in range(len(Mat1List)):
        ## count tells us how many of MatList[i] are in Mat1
        #print(str(i))
        Mat1=Mat1List[i]
        Mat2=Mat2List[i]
        
        if Mat1List[i].shape[0]==0:
            Mat1Str=" "
            subCount[i]=Mat2[~Mat2['proteinacc_start_stop_pre_post_.'].str.contains(Mat1Str)].shape[0]
            subMats[i]=Mat2[Mat2['proteinacc_start_stop_pre_post_.'].str.contains(Mat1Str)]
            print("Pep:"+str(len(Mat2[~Mat2['proteinacc_start_stop_pre_post_.'].str.contains(Mat1Str)]['Sequence'].unique())))
            print("UnqPep:"+str(len(temp[~temp['proteinacc_start_stop_pre_post_.'].str.contains(";")]['Sequence'].unique())))
        else:
            Mat1L=Mat1[column1].tolist()
            #colnames=Mat2.columns.values
            #print(colnames)
            idx=[]
            count=0
            Mat1Str="|".join(Mat1L)
            subCount[i]=Mat2[~Mat2['proteinacc_start_stop_pre_post_.'].str.contains(Mat1Str)].shape[0]
            print("Pep:"+str(len(Mat2[~Mat2['proteinacc_start_stop_pre_post_.'].str.contains(Mat1Str)]['Sequence'].unique())))
            temp=Mat2[~Mat2['proteinacc_start_stop_pre_post_.'].str.contains(Mat1Str)]
            print("UnqPep:"+str(len(temp[~temp['proteinacc_start_stop_pre_post_.'].str.contains(";")]['Sequence'].unique())))
            subMats[i]=Mat2[Mat2['proteinacc_start_stop_pre_post_.'].str.contains(Mat1Str)]
    return {'Count':subCount,'Matrix':subMats}
def mergedFile(res1,res1PSM, nonCDSPrtOutFileNameList,nonCDSPSMOutFileNameList):
    print("DirName:"+str(os.path.dirname(nonCDSPrtOutFileNameList[0])))
    mergedPrtFile=open(os.path.dirname(nonCDSPrtOutFileNameList[0])+"/nonCodingPrt.tsv",'w')
    mergedPSMFile=open(os.path.dirname(nonCDSPSMOutFileNameList[0])+"/nonCodingPSM.tsv",'w')
    for i in range(len(nonCDSPrtOutFileNameList)):
        res1['Matrix'][i].to_csv(nonCDSPrtOutFileNameList[i],sep='\t',index=False)
        res1PSM['Matrix'][i].to_csv(nonCDSPSMOutFileNameList[i],sep='\t',index=False)
        classCol=[]
        if i==0:
            classCol=['Non-sence mediated decay']
        elif i==1:
            classCol=['Processed transcript']
        elif i==2:
            classCol=['Retained intron']
        elif i==3:
            classCol=['Anti-sense']
        elif i==4:
            classCol=['lincRNA']
        prot=res1['Matrix'][i]
        psm=res1PSM['Matrix'][i]
        infoProt=pd.DataFrame(classCol*prot.shape[0],columns=['NonProteinCodingClass'],index=prot.index)
        infoPSM=pd.DataFrame(classCol*psm.shape[0],columns=['NonProteinCodingClass'],index=psm.index)
        prot=prot.join(infoProt)
        psm=psm.join(infoPSM)
        if i==0:
            prot.to_csv(mergedPrtFile,sep='\t',index=False)
            psm.to_csv(mergedPSMFile,sep='\t',index=False)
        else:
            prot.to_csv(mergedPrtFile,sep='\t', header=False, index=False)
            psm.to_csv(mergedPSMFile,sep='\t', header=False, index=False)
    mergedPSMFile.close()
    mergedPrtFile.close()

    
def main(knownPro, knownProSAP, iso, isoSAP, nonCDSFilenameList, nonCDSPSMFileNameList, nonCDSPrtOutFileNameList,nonCDSPSMOutFileNameList):
    ## Read a) known proteins, b) known proteins with SAPs including ALT/INDELs, c) isoforms, and d) isoforms with SAPs/ALT/INDELs
    
    ## Known Protein
    knownProteins=replaceCharacter(readList(knownPro,','),"ORF Id",';',',')
    ## Known with SAPs/ALT/INDELs
    knownProteinsSAPs=replaceCharacter(readList(knownProSAP,','),"ORF Id",';',',')
    ## Isoforms
    isoforms=replaceCharacter(readList(iso,','),"ORF Id",';',',')
    ## Isoforms with SAPs/ALT/INDELs
    isoformsSAPs=replaceCharacter(readList(isoSAP,','),"ORF Id",';',',')
    
    ## Read i) anisense, ii) linc, iii) nonsense-mediated decay, iv) processed transcript, v) reatined intron and any other ORFs from
    ## non CDS RNA bio type.
    nonCDS=[None]*len(nonCDSFilenameList)
    nonCDSPSMs=[None]*len(nonCDSPSMFileNameList)
    for i in range(len(nonCDSFilenameList)):
        nonCDS[i]=readList(nonCDSFilenameList[i],'\t')
        #print(nonCDS[i].columns.values)
        print("nonCDS:"+str(len(nonCDS[i])))
    
    for i in range(len(nonCDSPSMFileNameList)):
        nonCDSPSMs[i]=readList(nonCDSPSMFileNameList[i],'\t')
        #print(nonCDS[i].columns.values)
        print("nonCDSPSMs:"+str(len(nonCDSPSMs[i])))
        print("nonCDSPep:"+str(len(nonCDSPSMs[i]['Sequence'].unique())))
        print("nonCDSUniqPep:"+str(len(nonCDSPSMs[i][~nonCDSPSMs[i]['proteinacc_start_stop_pre_post_.'].str.contains(";")]['Sequence'].unique())))

    ##MSGF+ protein file has comma in protein accession, whereas MatList files, i.e. known proteins etc. have ';'.
    ## Remove all proteins/ORFs from i-v if they are found in a. Make separate lists of those ORFs/proteins from the i-v lists that
    ## exist in b, c or d.
    print("Known Protein")
    res1=subsetProteins(knownProteins,nonCDS,'ORF Id')
    print("res1"+str(res1['Count']))
    res1PSM=subsetPSMs(res1['Matrix'], nonCDSPSMs, 'protein.accession')
    print("Known Protein with Var")
    mergedFile(res1,res1PSM, nonCDSPrtOutFileNameList,nonCDSPSMOutFileNameList)
   
    res2=subsetProteins(knownProteinsSAPs,res1['Matrix'],'ORF Id')
    print("res2"+str(res2['Count']))
    res2PSM=subsetPSMs(res2['Matrix'], res1PSM['Matrix'], 'protein.accession')
    
    print("Isoform")
    res3=subsetProteins(isoforms,res2['Matrix'],'ORF Id')
    print("res3"+str(res3['Count']))
    res3PSM=subsetPSMs(res3['Matrix'], res2PSM['Matrix'], 'protein.accession')
    
    print("Isoform with Var")
    res4=subsetProteins(isoformsSAPs,res3['Matrix'],'ORF Id')
    print("res4"+str(res4['Count']))
    #print("res4 Antisense:")
    #print(res4['Matrix'][3]['description'])
    res4PSM=subsetPSMs(res4['Matrix'], res3PSM['Matrix'], 'protein.accession')
    
    print("res1PSM:"+str(res1PSM['Count']))
    print("res2PSM:"+str(res2PSM['Count']))
    print("res3PSM:"+str(res3PSM['Count']))
    print("res4PSM:"+str(res4PSM['Count']))
    
    print("NMD,PT,RI,ANTI,LINC")
    st1=""
    st2=""
    st3=""
    st4=""
    for i in range(len(res1PSM['Matrix'])):
        st1=st1+",Pep:"+str(len(res1PSM['Matrix'][i]['Sequence'].unique()))
        st2=st2+",Pep:"+str(len(res2PSM['Matrix'][i]['Sequence'].unique()))
        st3=st3+",Pep:"+str(len(res3PSM['Matrix'][i]['Sequence'].unique()))
        st4=st4+",Pep:"+str(len(res4PSM['Matrix'][i]['Sequence'].unique()))
        st1=st1+"|UniqPep:"+str(len(res1PSM['Matrix'][i][~res1PSM['Matrix'][i]['proteinacc_start_stop_pre_post_.'].str.contains(";")]["Sequence"].unique()))
        st2=st2+"|UniqPep:"+str(len(res2PSM['Matrix'][i][~res2PSM['Matrix'][i]['proteinacc_start_stop_pre_post_.'].str.contains(";")]["Sequence"].unique()))
        st3=st3+"|UniqPep:"+str(len(res3PSM['Matrix'][i][~res3PSM['Matrix'][i]['proteinacc_start_stop_pre_post_.'].str.contains(";")]["Sequence"].unique()))
        st4=st4+"|UniqPep:"+str(len(res4PSM['Matrix'][i][~res4PSM['Matrix'][i]['proteinacc_start_stop_pre_post_.'].str.contains(";")]["Sequence"].unique()))
    print(st1)
    print(st2)
    print(st3)
    print(st4)
    
    return [res1,res2,res3,res4,res1PSM,res2PSM,res3PSM,res4PSM]


knownPro="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsV7.csv"
knownProSAP="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsSAPsV7.csv"
iso="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsV7.csv"
isoSAP="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsSAPsV7.csv"

nonCDSPrtFileNameList=["D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonsense_mediated_decayPrt.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/processed_transcriptPrt.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/retained_intronPrt.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/antisensePrt.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/lincPrt.tsv"]

nonCDSPrtOutFileNameList=["D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonsense_mediated_decayPrt2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/processed_transcriptPrt2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/retained_intronPrt2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/antisensePrt2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/lincPrt2.tsv"]

nonCDSPSMFileNameList=["D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonsense_mediated_decay.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/processed_transcript.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/retained_intron.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/antisense.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/linc.tsv"]
nonCDSPSMOutFileNameList=["D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonsense_mediated_decay2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/processed_transcript2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/retained_intron2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/antisense2.tsv",
                       "D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/linc2.tsv"]
nonCDSFilenameList=nonCDSPrtFileNameList
allres=main(knownPro, knownProSAP, iso, isoSAP, nonCDSPrtFileNameList,nonCDSPSMFileNameList,nonCDSPrtOutFileNameList,nonCDSPSMOutFileNameList)
'''
print("All Res:"+str(allres[3]['Count']))
for i in range(len(allres[3]['Matrix'])):
    print("last Mat dimensions:"+str(allres[0]['Matrix'][i].shape))
    print("last Mat dimensions:"+str(allres[4]['Matrix'][i].shape))
    #print("Ind Pep Count:"+str(allres[4]['Count']))
    print("Pep Count:"+str(len(allres[6]['Matrix'][i]['Sequence'].unique())))

for i in range(int(len(allres)/2)):
    print("Protein class:"+str(i))
    for j in range(1,len(allres[i]['Matrix'])):
        print("\tnonCDS Class:"+str(j))
        if j<int(len(allres)/2)-1:
            print("\t\tProtein Count:"+str(allres[i]['Matrix'][j].shape[0]-allres[i]['Matrix'][j+1].shape[0]))
        else:
            print("\t\tProtein Count:"+str(allres[i]['Matrix'][j].shape[0]))
'''
'''
for i in range(int(len(allres)/2),len(allres)):
    print("Protein class:"+str(i))
    for j in range(1,len(allres[i]['Matrix'])):
        print("\tnonCDS Class:"+str(j))
        if j<len(allres[i]['Matrix'])-1:
            protAccJ=allres[i]['Matrix'][j]
            protAccJ1=allres[i]['Matrix'][j+1]
            uniqJ=len(protAccJ[~protAccJ['proteinacc_start_stop_pre_post_.'].str.contains(";")]["Sequence"].unique())
            uniqJ1=len(protAccJ1[~protAccJ1['proteinacc_start_stop_pre_post_.'].str.contains(";")]["Sequence"].unique())
            print("\t\tPeptide Count:"+str((len(allres[i]['Matrix'][j]['Sequence'].unique()))-(len(allres[i]['Matrix'][j+1]['Sequence'].unique()))))
            print("\t\tPSMs Count:"+str(allres[i]['Matrix'][j].shape[0]-allres[i]['Matrix'][j+1].shape[0]))
            print("\t\tUnique Peptide Count:"+str(uniqJ-uniqJ1))
            print("\t\tPep Count Inc:"+str(len(allres[i]['Matrix'][j]['Sequence'].unique())))
        else:
            protAccJ=allres[i]['Matrix'][j]
            uniqJ=len(protAccJ[~protAccJ['proteinacc_start_stop_pre_post_.'].str.contains(";")]["Sequence"].unique())
            print("\t\tPeptide Count:"+str(len(allres[i]['Matrix'][j]['Sequence'].unique())))
            print("\t\tPSMs Count:"+str(allres[i]['Matrix'][j].shape[0]))
            print("\t\tUnique Peptide Count:"+str(uniqJ))
            print("\t\tPep Count Inc:"+str(len(allres[i]['Matrix'][j]['Sequence'].unique())))
'''

'''
    NMD,PT,RI,ANTI,LINC
KP[10, 2, 5, 0, 0]
KPV[2, 0, 4, 0, 0]
ISO+V[24, 11, 22, 0, 2]
UNK[3, 5, 4, 1, 2]
res1PSM:[562, 327, 387, 6, 9]
res2PSM:[522, 327, 331, 6, 9]
res3PSM:[181, 81, 31, 6, 6]
res4PSM:[6, 17, 5, 6, 6]

nonCDS:40
nonCDS:20
nonCDS:36
nonCDS:1
nonCDS:4
nonCDSPSMs:692
nonCDSPep:194
nonCDSUniqPep:82
nonCDSPSMs:376
nonCDSPep:80
nonCDSUniqPep:22
nonCDSPSMs:402
nonCDSPep:153
nonCDSUniqPep:19
nonCDSPSMs:6
nonCDSPep:2
nonCDSUniqPep:0
nonCDSPSMs:9
nonCDSPep:4
nonCDSUniqPep:0
NMD,PT,RI,ANTI,LINC
[1]KP+[2,3,4]   [Pep:141|UniqPep:51,Pep:51|UniqPep:0,Pep:143|UniqPep:15,Pep:2|UniqPep:0,Pep:4|UniqPep:0
[2]KPV+[3,4]  [Pep:137|UniqPep:51,Pep:51|UniqPep:0,Pep:124|UniqPep:15,Pep:2|UniqPep:0,Pep:4|UniqPep:0
[3]ISO+V+[4][Pep:12|UniqPep:0,Pep:21|UniqPep:0,Pep:14|UniqPep:0,Pep:2|UniqPep:0,Pep:2|UniqPep:0
[4]UNK  [Pep:2|UniqPep:0,Pep:5|UniqPep:0,Pep:2|UniqPep:0,Pep:2|UniqPep:0,Pep:2|UniqPep:0


'''
### This Python code reads variations of identified ORFs (reverse and peptide 2 threshold applied) and variations with
### peptide evidences. It then counts number of ORFs in there and count number of variation for the analysis
### purpose.
import pandas as pd

def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj
def nonProteinCoding(fileNames,sep):
    for i in range(len(fileNames)):
        ##Read files
        print(fileNames[i])

def readAndSplitInfo(vcfFileName,flag):
    vcf=readFile(vcfFileName,'\t')
    
    ## vcf file fields: Chr, POS within Protein,ID, REF, ALT, INFO(SubjectId=P09417-2;QueryId=Dataset_A_asmbl_41426_ORF20_Frame_3_84-446;Alignment=[QueryLength=121:QueryStart=1:QueryEnd=116:SubjectLength=213:SubjectStart=1:SubjectEnd=116];Type:SSAP;QPOS:116)
    ## split the INFO column and add them to main data frame
    if flag==1:
        info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS'])
    else:
        info=pd.DataFrame(vcf.INFO.str.split(';').tolist(),columns=['SubjectID','QueryID','Alignment','Type','QPOS','PeptideCount','UniquePeptideCount','Peptides','Score'])
    vcf=vcf.drop('INFO',1)
    vcf=vcf.join(info)
    vcf.SubjectID=vcf.SubjectID.str.replace('SubjectId=','')
    vcf.QueryID=vcf.QueryID.str.replace('QueryId=','')
    #print(vcf.QueryID[0:5])
    vcf.Alignment=vcf.Alignment.str.replace('Alignment=\[','')
    vcf.Alignment=vcf.Alignment.str.replace('\]','')
    vcf.Type=vcf.Type.str.replace('Type:','')
    vcf.QPOS=vcf.QPOS.str.replace('QPOS:','')
    return vcf
    
def subsettingIdentifiedORFVariation(vcfFileNameIdent, vcfFileNameEvd):
    ##read vcf file
    vcfIdent=readAndSplitInfo(vcfFileNameIdent,1)
    vcfEvd=readAndSplitInfo(vcfFileNameEvd,0)
    ##groups vcf entries by ORF/Query ids.
    vcfIdentGrouped=vcfIdent.groupby('QueryID')
    vcfEvdGrouped=vcfEvd.groupby('QueryID')
    ##Identified
    print(":::::Identified:::::")
    print("Total number of ORFs:"+str(len(vcfIdentGrouped)))
    print("Total numbber of SAAVs/SAPs [SSAPs+SAPs]:"+str(vcfIdent["Type"].str.contains("SAP|SSAP").sum())+" ["
          +str(vcfIdent["Type"].str.contains("SSAP").sum())+"+"
          +str(vcfIdent["Type"].str.contains("SAP").sum())+"]")
    print("Total number of ALTs [SALT/ALT]:"+str(vcfIdent["Type"].str.contains("ALT|SALT").sum())+" ["
          +str(vcfIdent["Type"].str.contains("SALT").sum())+"+"
          +str(vcfIdent["Type"].str.contains("ALT").sum())+"]")
    print("Total number of Insertion:"+str(vcfIdent["Type"].str.contains("INS").sum()))
    print("Total number of Deletion:"+str(vcfIdent["Type"].str.contains("DEL").sum()))
    
    ##Peptide Evidence
    print(":::::Peptide Evidence:::::")
    print("Total number of ORFs:"+str(len(vcfEvdGrouped)))
    print("Total numbber of SAAVs/SAPs [SSAPs+SAPs]:"+str(vcfEvd["Type"].str.contains("SAP|SSAP").sum())+" ["
          +str(vcfEvd["Type"].str.contains("SSAP").sum())+"+"
          +str(vcfEvd["Type"].str.contains("SAP").sum())+"]")
    print("Total number of ALTs [SALT/ALT]:"+str(vcfEvd["Type"].str.contains("ALT|SALT").sum())+" ["
          +str(vcfEvd["Type"].str.contains("SALT").sum())+"+"
          +str(vcfEvd["Type"].str.contains("ALT").sum())+"]")
    print("Total number of Insertion:"+str(vcfEvd["Type"].str.contains("INS").sum()))
    print("Total number of Deletion:"+str(vcfEvd["Type"].str.contains("DEL").sum()))

identFile="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IdentifiedORFsVariationPep2.vcf"
evdFile="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_with_Location_VariationV7PeptideEvidencePep2.vcf"

subsettingIdentifiedORFVariation(identFile,evdFile)
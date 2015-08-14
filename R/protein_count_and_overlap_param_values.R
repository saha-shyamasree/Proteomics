################ These are the parameters for protein_count_and_overlap.
f1,f2,f3,d1,d2,d3,rev,peptide,pepThreshold,upper, Mat1FileNamePostfix, Mat2FileNamePostfix

########################################## These are values ##############################################################################
#################### Common Values ###########
upper=0.000000000000000000000000000001
pepThreshold=1
rev=1
peptide=1
########################################## 06-08-2015 trinityV5 vs PASA #############################################################################

pPrtFile="pasa_assemblyV1+fdr+th+grouping+prt.csv"
tPrtFile="trinity_PITORF+fdr+th+grouping+prtV5.csv"
pPrtDir="D:/data/Results/Human-Adeno/Identification/PASA/sORF/"
tPrtDir="D:/data/Results/Human-Adeno/Identification/sORF/"
pBlastFile="human_adeno_mydb_pasa.assemblies_ORFsV1.csv"
tBlastFile="trinityV5Match.csv"
pblastDir="D:/data/blast/blastCSV/PASA/Human-Adeno/"
tblastDir="D:/data/blast/blastCSV/sORF/"
Mat1FileNamePostfix="TrinityV5"
Mat2FileNamePostfix="PASA"
PAG_overlap_ThirdDB_Map_Main(tPrtFile,pPrtFile,tBlastFile,pBlastFile,tPrtDir,pPrtDir,tblastDir,pblastDir,rev,peptide,pepThreshold,upper, Mat1FileNamePostfix, Mat2FileNamePostfix)


########################################## 04-08-2015, trinityV5 uniprot overlap #########################################################

f1="trinity_PITORF+fdr+th+grouping+prtV5.csv"
f2="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
f3="trinityV5Match.csv"
d1="D:/data/Results/Human-Adeno/Identification/sORF/" 
d2="D:/data/Results/Human-Adeno/Identification/"
#d1=d2
d3="D:/data/blast/blastCSV/sORF/"
Mat1FileNamePostfix="TrinityV5"
Mat2FileNamePostfix="Uniprot"
PAG_overlap_Main(f1,f2,f3,d1,d2,d3,rev,peptide,pepThreshold,upper, Mat1FileNamePostfix, Mat2FileNamePostfix)

########################################## 18-05-2015, PASA uniprot overlap #########################################################
f1="pasa_assemblyV1+fdr+th+grouping+prt.csv"
f2="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
f3="human_adeno_mydb_pasa.assemblies_ORFsV1.csv"
d1="D:/data/Results/Human-Adeno/Identification/PASA/sORF/"
d2="D:/data/Results/Human-Adeno/Identification/"
#d1=d2
d3="D:/data/blast/blastCSV/PASA/Human-Adeno/"
Mat1FileNamePostfix="PASA"
Mat2FileNamePostfix="Uniprot"
upper=0.000000000000000000000000000001
PAG_overlap_Main(f1,f2,f3,d1,d2,d3,rev,peptide,pepThreshold,upper, Mat1FileNamePostfix, Mat2FileNamePostfix)

#####################################################################################################################################
filenameBHT="trinity_bat+fdr+th+grouping+prt.csv"
filenameBHU="uniprot_bat+fdr+th+grouping+prt.csv"
folder="C:/Users/shyama/Dropbox/results/Bat_human_hendra/Bat/"
TBH=proteinGroupFiltered(proteinGroup(filenameBHT,folder),rev,peptide,pepThreshold)
UBH=proteinGroupFiltered(proteinGroup(filenameBHU,folder),rev,peptide,pepThreshold)
filenameDB="uni_db_p_alecto_blast.xml2.csv"
dbfolder="C:/Users/shyama/Dropbox/results/blastdb/blastCSV/"
blastBH=blastFilterEval(blastFilterMatch(blast(filenameDB,dbfolder),1),upper)
TBH=replaceComma(TBH,3)
TBH=upto(TBH,3,',')
blastBH=upto(blastBH,4,' ')
UBH[,'protein.accession']=sub("^sw","sp",UBH[,'protein.accession'])
TBH[,'protein.accession']=replaceIds(TBH,blastBH)

length(overlap(proteinGrpAnchor(TBH)[,'protein.accession'],UBH))
length(grep("_HENDH",unique(proteinGrpAnchor(TBH)[,'protein.accession'])))

BTExc=proteinGrpAnchor(TBH)[exclusive(proteinGrpAnchor(TBH)[,'protein.accession'],UBH),c('protein.accession','distinct.peptide.sequences')]
length(which(as.numeric(BTExc[grep("^comp",BTExc[,'protein.accession']),'distinct.peptide.sequences'])>2))
length(which(as.numeric(BTExc[,'distinct.peptide.sequences'])>2))

HUH_hendra=proteinGrpAnchor(HUHendra)[grep("_HENDH",proteinGrpAnchor(HUHendra)[,'protein.accession']),'protein.accession']
BUH_hendra=proteinGrpAnchor(UBH)[grep("_HENDH",proteinGrpAnchor(UBH)[,'protein.accession']),'protein.accession']

######################## Trinity vs Uniprot  #################################
twoDataSetsOverlap(f1,f2,f3,d1,d2,d3,1,1,1,0.000000000000000000000000000001)

f1="trinity_PITORF+fdr+th+grouping+prt.csv"
f2="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
f3="trinity_PITORF_human_adeno_blast2.csv"
#d1="C:/Users/shyama/Dropbox/results/Human_adenovirus/"
d2="D:/data/Results/Human-Adeno/Identification/"
d1=d2
d3="D:/data/blast/blastCSV/"
Mat1FileNamePostfix="Trinity"
Mat2FileNamePostfix="Uniprot"
PAG_overlap_Main(f1,f2,f3,d1,d2,d3,rev,peptide,pepThreshold,upper, Mat1FileNamePostfix, Mat2FileNamePostfix)

venn.plot <- venn.diagram(list(Uniprot = 1:3015, ORFs = 247:3173),"Venn_2set_ProteinTrinity.tiff",fill=c("black","white"), alpha=c(0.65,.75), offset=0.5,euler.d=TRUE,scaled=TRUE, ext.text=TRUE, cex=c(2,2,2), cat.col=c("white","white"), ext.line.lwd=0, label.col=c("white","black","white"), ext.percent=c(.25,.25,.25))

######################### Cufflink with uniprot inc Adeno Comparison  ############################################

#f1="cufflinks_PITORF+fdr+th+grouping+prt.csv"
f1="cufflinks_main_PITORF+fdr+th+grouping+prt.csv"
f2="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
#f3="cufflinks-ORF.csv"
f3="cufflinks_main-ORF.csv"
#d1="C:/Users/shyama/Dropbox/results/Human_adenovirus/"
d2="D:/data/Results/Human-Adeno/Identification/"
d1=d2
d3="D:/data/blast/blastCSV/"

venn.plot <- venn.diagram(list(Uniprot = 1:2986, ORFs = 551:3264),"Venn_2set_ProteinCuff.tiff",fill=c("black","white"), alpha=c(0.65,.75), offset=0.5,euler.d=TRUE,scaled=TRUE, ext.text=TRUE, cex=c(2,2,2), cat.col=c("white","white"), ext.line.lwd=0, label.col=c("white","black","white"), ext.percent=c(.25,.25,.25))

######################### Cufflink with uniprot human Comparison  ############################################

f1="cufflinks_main_PITORF+fdr+th+grouping+prt.csv"
f2="human_uniprot+fdr+th+grouping+prt.csv"
f3="human_cufflinks_mainORF.csv"
#d1="C:/Users/shyama/Dropbox/results/Human_adenovirus/"
d2="D:/data/Results/Human-Adeno/Identification/"
d1=d2
d3="D:/data/blast/blastCSV/"


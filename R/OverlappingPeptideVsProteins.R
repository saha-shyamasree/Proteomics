##this code compares the PIT search results against the Standard search. Aim of this code is to check how many search DB specific peptide identification
##contribute towards search DB specific protein identification.

source("D:/Code/Proteomics/R/RLib.R")



dir1="D:/data/Results/Human-Adeno/GIOPaperResults/trinity/"
dir2="D:/data/Results/Human-Adeno/GIOPaperResults/Uniprot-trinity/"
dir3="D:/data/Results/Human-Adeno/GIOPaperResults/"
dir4=dir3
fname1="trinity_PITORF+fdr+th+grouping.csv"
fname1prt="trinity_PITORF+fdr+th+grouping+prt.csv"

fname2="DM_from_raw_uniprot+fdr+th+grouping.csv"
fname2prt="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
fname3="PAGComparison_TrinityOnly_withUniprot.tsv"
fname4="AnchortoAnchor_PAGs_uni_trinity1e-30.tsv"
fname5="AnchortoSub_PAGs_uni_trinity1e-30.tsv"

##peptide identification
Mat1PSM=readPSM(fname1,dir1)
Mat1Prt=proteinGroup(fname1prt,dir1)


##Peptide only identified for PIT-without genome search and their association with one hit wonder proteins
Mat1PrtSingle=oneHitWonderProteinGrp(Mat1Prt,1)
Mat1PrtSinglePep=c(Mat1PrtSingle[,'unique.peptides'],Mat1PrtSingle[,'razor.peptides'])
Mat1PrtSinglePepStr=gsub("\\s+","",gsub(";+",';',paste(Mat1PrtSinglePep,collapse=";")))
Mat1PrtSinglePepUniq=unique(strsplit(Mat1PrtSinglePepStr,';')[[1]])
head(Mat1PrtSinglePepUniq)
##if "" is there then do following
Mat1PrtSinglePepUniq=Mat1PrtSinglePepUniq[-1]

##How many of these peptides are from non-overlapping peptides?
length(which(unique(Mat1PepUniq[,'Sequence']) %in% Mat1PrtSinglePepUniq))


##Peptides supporting PIT proteins only.
Mat1Exc=proteinList(fname3,dir3,"\t")
Mat1ExcPep=c(Mat1Exc[,'unique.peptides'],Mat1Exc[,'razor.peptides'])
Mat1ExcPepStr=gsub("\\s+","",gsub(";+",';',paste(Mat1ExcPep,collapse=";")))
Mat1ExcPepUniq=unique(strsplit(Mat1ExcPepStr,';')[[1]])

Mat1AtoA=proteinList(fname4,dir3,"\t")
Mat1AtoAPep=c(Mat1AtoA[,'unique.peptides'],Mat1AtoA[,'razor.peptides'])
Mat1AtoAPepStr=gsub("\\s+","",gsub(";+",';',paste(Mat1AtoAPep,collapse=";")))
Mat1AtoAPepUniq=unique(strsplit(Mat1AtoAPepStr,';')[[1]])

Mat1AtoSub=proteinList(fname5,dir3,"\t")
Mat1AtoSubPep=c(Mat1AtoSub[,'unique.peptides'],Mat1AtoSub[,'razor.peptides'])
Mat1AtoSubPepStr=gsub("\\s+","",gsub(";+",';',paste(Mat1AtoSubPep,collapse=";")))
Mat1AtoSubPepUniq=unique(strsplit(Mat1AtoSubPepStr,';')[[1]])

fname6="PAGComparison_UniprotOnly_withTrinity.tsv"


Mat2PSM=readPSM(fname2,dir2)
Mat2Prt=proteinGroup(fname2prt,dir2)

##Peptide only identified for Standard search and their association with one hit wonder proteins
Mat2PrtSingle=oneHitWonderProteinGrp(Mat2Prt,1)
Mat2PrtSinglePep=c(Mat2PrtSingle[,'unique.peptides'],Mat2PrtSingle[,'razor.peptides'])
Mat2PrtSinglePepStr=gsub("\\s+","",gsub(";+",';',paste(Mat2PrtSinglePep,collapse=";")))
Mat2PrtSinglePepUniq=unique(strsplit(Mat2PrtSinglePepStr,';')[[1]])
head(Mat2PrtSinglePepUniq)
##if "" is there then do following
Mat2PrtSinglePepUniq=Mat2PrtSinglePepUniq[-1]

##How many of these peptides are from non-overlapping peptides?
length(which(unique(Mat2PepUniq[,'Sequence']) %in% Mat2PrtSinglePepUniq))


##Peptides supporting Standard only proteins
Mat2Exc=proteinList(fname6,dir3,"\t")
Mat2ExcPep=c(Mat2Exc[,'unique.peptides'],Mat1Exc[,'razor.peptides'])
Mat2ExcPepStr=gsub("\\s+","",gsub(";+",';',paste(Mat2ExcPep,collapse=";")))
Mat2ExcPepUniq=unique(strsplit(Mat2ExcPepStr,';')[[1]])


column1=11
column2=2



##Peptides only identified for PIT search
Mat1PepUniq=Mat1PSM[-overlapPepIndices(Mat1PSM, Mat2PSM, column1, column1),]
Mat1PepUniqSeq=unique(Mat1PepUniq[,'Sequence'])

##Overlaping peptides.
Mat1PepOver=Mat1PSM[overlapPepIndices(Mat1PSM, Mat2PSM, column1, conamelumn1),]

###Peptides only identified for PIT search association with PIT only proteins
length(which(unique(Mat1ExcPepUniq) %in% unique(Mat1PepUniq)))

###Peptides only identified for PIT search, association with overlaping Anchor proteins' peptides.
length(unique(Mat1PepUniq[which(Mat1PepUniq[,'Sequence'] %in% Mat1AtoAPepUniq),'Sequence']))
length(which(unique(Mat1PepUniq[,'Sequence']) %in% Mat1AtoAPepUniq))

##Peptides only identified for Standard search
Mat2PepUniq=Mat2PSM[-overlapPepIndices(Mat2PSM, Mat1PSM, column1, column1),]
Mat2PepUniqSeq=unique(Mat2PepUniq[,'Sequence'])

##Overlaping peptides.
Mat2PepOver=Mat2PSM[overlapPepIndices(Mat2PSM, Mat1PSM, column1, column1),]

###Peptides only identified for Standard search association with PIT only proteins
length(which(unique(Mat2ExcPepUniq) %in% unique(Mat2PepUniq)))

###Peptides only identified for Standard search, association with overlaping Anchor proteins' peptides.
length(unique(Mat2PepUniq[which(Mat2PepUniq[,'Sequence'] %in% Mat1AtoAPepUniq),'Sequence']))


###################### Human adeno virus Standard vs PIT genome guided  ##########################

dirc="D:/data/Results/Human-Adeno/GIOPaperResults/cufflinks/"
diru="D:/data/Results/Human-Adeno/GIOPaperResults/Uniprot-Cufflinks/"
dircomm="D:/data/Data/fasta/GIOPaper/"

fnamec="cufflinks_main_PITORF+fdr+th+grouping.csv"
fnamecprt="cufflinks_main_PITORF+fdr+th+grouping+prt.csv"
fnameu="human_uniprot+fdr+th+grouping.csv"
fnameuprt="human_uniprot+fdr+th+grouping+prt.csv"
fnamec3="PAGComparison_CuffOnly_withUniprot1e-30.tsv"
fnamec4="AnchortoAnchor_PAGs_uni_cuff1e-30.tsv"
fnamec5="AnchortoSub_PAGs_uni_cuff1e-30.tsv"

fnameu3="PAGComparison_UniprotOnly_withCuff1e-30.tsv"
##following 2 need correction
fnameu4="AnchortoAnchor_PAGs_uni_cuff1e-30.tsv"
fnameu5="AnchortoSub_PAGs_uni_cuff1e-30.tsv"

column1=11
column2=2

MatCPSM=readPSM(fnamec,dirc)
MatCPrt=proteinGroup(fnamecprt,dirc)
##Peptides supporting PIT-genome guided proteins only.
MatCExc=proteinList(fnamec3,dir3,"\t")
MatCExcPep=c(MatCExc[,'unique.peptides'],MatCExc[,'razor.peptides'])
MatCExcPepStr=gsub("\\s+","",gsub(";+",';',paste(MatCExcPep,collapse=";")))
MatCExcPepUniq=unique(strsplit(MatCExcPepStr,';')[[1]])

##Peptides supporting Standard-genome guided proteins only.
MatUExc=proteinList(fnameu3,dir3,"\t")
MatUExcPep=c(MatUExc[,'unique.peptides'],MatUExc[,'razor.peptides'])
MatUExcPepStr=gsub("\\s+","",gsub(";+",';',paste(MatUExcPep,collapse=";")))
MatUExcPepUniq=unique(strsplit(MatUExcPepStr,';')[[1]])


##Peptides only identified for PIT-genome guided search
MatCPepUniq=MatCPSM[-overlapPepIndices(MatCPSM, MatUPSM, column1, column1),]
##Overlaping peptides.
MatCPepOver=MatCPSM[overlapPepIndices(MatCPSM, MatUPSM, column1, column1),]

###Peptides only identified for PIT-genome guided search association with PIT-genome guided only proteins
length(which(unique(MatCExcPepUniq) %in% unique(MatCPepUniq[,'Sequence'])))

###Peptides only identified for PIT search, association with overlaping Anchor proteins' peptides.
length(unique(MatCPepUniq[which(MatCPepUniq[,'Sequence'] %in% MatCAtoAPepUniq),'Sequence']))

##Peptide only identified for PIT-genome guided search and their association with one hit wonder proteins
MatCPrtSingle=oneHitWonderProteinGrp(MatCPrt,1)
MatCPrtSinglePep=c(MatCPrtSingle[,'unique.peptides'],MatCPrtSingle[,'razor.peptides'])
MatCPrtSinglePepStr=gsub("\\s+","",gsub(";+",';',paste(MatCPrtSinglePep,collapse=";")))
MatCPrtSinglePepUniq=unique(strsplit(MatCPrtSinglePepStr,';')[[1]])
head(MatCPrtSinglePepUniq)
##if "" is there then do following
MatCPrtSinglePepUniq=MatCPrtSinglePepUniq[-1]
MatCPepUniqSeq=unique(MatCPepUniq[,'Sequence'])
##How many of these peptides are from non-overlapping peptides?
length(which( MatCPepUniqSeq %in% MatCPrtSinglePepUniq))
##peptides that are contributing to one hit wonders.
MatCPepUniqSeqSingle=MatCPepUniqSeq[which(MatCPepUniqSeq %in% MatCPrtSinglePepUniq)]

##Peptides only identified for Standard search
MatUPepUniq=MatUPSM[-overlapPepIndices(MatUPSM, MatCPSM, column1, column1),]
##Overlaping peptides.
MatUPepOver=MatUPSM[overlapPepIndices(MatUPSM, MatCPSM, column1, column1),]

###Peptides only identified for Standard search association with PIT only proteins
length(which(unique(MatUExcPepUniq) %in% unique(MatUPepUniq)))

MatUPSM=readPSM(fnameu,diru)
MatUPrt=proteinGroup(fnameuprt,diru)

##Peptide only identified for Standard-genome guided search and their association with one hit wonder proteins
MatUPrtSingle=oneHitWonderProteinGrp(MatUPrt,1)
MatUPrtSinglePep=c(MatUPrtSingle[,'unique.peptides'],MatUPrtSingle[,'razor.peptides'])
MatUPrtSinglePepStr=gsub("\\s+","",gsub(";+",';',paste(MatUPrtSinglePep,collapse=";")))
MatUPrtSinglePepUniq=unique(strsplit(MatUPrtSinglePepStr,';')[[1]])
head(MatUPrtSinglePepUniq)
##if "" is there then do following
MatUPrtSinglePepUniq=MatUPrtSinglePepUniq[-1]

MatUPepUniqSeq=unique(MatUPepUniq[,'Sequence'])
##How many of these peptides are from non-overlapping peptides?
length(which(MatUPepUniqSeq %in% MatUPrtSinglePepUniq))
MatUPepUniqSeqSingle=MatUPepUniqSeq[which(MatUPepUniqSeq %in% MatUPrtSinglePepUniq)]
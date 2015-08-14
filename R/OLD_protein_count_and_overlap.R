############ Old protein overlap funtions.

twoDataSetsExactOverlap <- function(f1,f2,f3,d1,d2,d3,rev=1,peptide=1,pepThreshold=1,gthreshold=1,lthreshold=1,rthreshold=1)
{
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    Mat2=proteinGroupFiltered(proteinGroup(f2,d1),rev,peptide,pepThreshold)
    blastDB=blastFilterLengthRatioMatch(blastFilterLongMatch(blastFilterGoodMatch(blastFilterMatch(blast(f3,d3),1),gthreshold),lthreshold),rthreshold)
    print("Mat1")
    print(dim(Mat1))
    Mat1=replaceComma(Mat1,3)
    print(dim(Mat1))
    Mat1=upto(Mat1,3,';')
    print(dim(Mat1))
    Mat1=upto(Mat1,3,' ')
    print(dim(Mat1))
    Mat2=replaceComma(Mat2,3)
    Mat2=upto(Mat2,3,';')
    Mat2=upto(Mat2,3,' ')
    blastDB=upto(blastDB,1,';')
    blastDB=upto(blastDB,1,',')
    blastDB=upto(blastDB,4,' ')
    Mat1[,'protein.accession']=sub("^sw","sp",Mat1[,'protein.accession'])
    print(dim(Mat1))
    Mat2[,'protein.accession']=sub("^sw","sp",Mat2[,'protein.accession'])
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    Mat2[,'protein.accession']=replaceIds(Mat2,blastDB)


    #print(dim(Mat1))
    #print("Mat1 overlap dim")
    #print(dim(Mat1))
    #write.table(Mat1,file=paste(d1,"tempMat1.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat2,file=paste(d2,"tempMat2.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #print(length(overlapIndices(proteinGrpAnchor(Mat1)[,'protein.accession'],Mat2)))
    Mat1Anchor=proteinGrpAnchor(Mat1)
    Mat2Anchor=proteinGrpAnchor(Mat2)
    print(paste("Mat1 Anchor protein:",length(unique(Mat1Anchor[,'protein.accession']))))
    print(paste("Mat2 Anchor protein:",length(unique(Mat2Anchor[,'protein.accession']))))
    print(paste("Overlap:Mat1 to Mat2",length(overlap(Mat1Anchor[,'protein.accession'],Mat2))))
    print(paste("Overlap:Mat2 to Mat1",length(overlap(Mat2Anchor[,'protein.accession'],Mat1))))
    print(paste("Overlap:Anchor to Anchor",length(overlap(Mat2Anchor[,'protein.accession'],Mat1Anchor))))
    #write.table(Mat1Anchor[overlapIndices(Mat1Anchor[,'protein.accession'],Mat2),],file=paste(d1,"overlap_",gthreshold,"_",lthreshold,"_",rthreshold,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat1Anchor[exclusive(Mat1Anchor[,'protein.accession'],Mat2),],file=paste(d1,"only",f1,"_",gthreshold,"_",lthreshold,"_",rthreshold,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat2Anchor[exclusive(Mat2Anchor[,'protein.accession'],Mat1),],file=paste(d2,"only",f2,"_",gthreshold,"_",lthreshold,"_",rthreshold,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #list(Mat1_back,Mat2_back)
}

twoDataSetsOverlap <- function(f1,f2,f3,d1,d2,d3,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    Mat2=proteinGroupFiltered(proteinGroup(f2,d2),rev,peptide,pepThreshold)
    blastDB=blastFilterEval(blastFilterMatch(blast(f3,d3),1),upper)
    print("Mat1")
    print(dim(Mat1))
    Mat1=replaceComma(Mat1,3)
    print(dim(Mat1))
    Mat1=upto(Mat1,3,';')
    print(dim(Mat1))
    Mat1=upto(Mat1,3,' ')
    print(dim(Mat1))
    Mat2=replaceComma(Mat2,3)
    Mat2=upto(Mat2,3,';')
    Mat2=upto(Mat2,3,' ')
    blastDB=upto(blastDB,1,';')
    blastDB=upto(blastDB,1,',')
    blastDB=upto(blastDB,4,' ')
    Mat1[,'protein.accession']=sub("^sw","sp",Mat1[,'protein.accession'])
    print(dim(Mat1))
    Mat2[,'protein.accession']=sub("^sw","sp",Mat2[,'protein.accession'])
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    Mat2[,'protein.accession']=replaceIds(Mat2,blastDB)

    #print(dim(Mat1))
    #print("Mat1 overlap dim")
    #print(dim(Mat1))
    #write.table(Mat1,file=paste(d1,"tempMat1.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat2,file=paste(d2,"tempMat2.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #print(length(overlapIndices(proteinGrpAnchor(Mat1)[,'protein.accession'],Mat2)))
    Mat1Anchor=proteinGrpAnchor(Mat1)
    Mat2Anchor=proteinGrpAnchor(Mat2)
    Mat2Area=length(unique(Mat2Anchor[,'protein.accession']))
    Mat1Area=length(unique(Mat1Anchor[,'protein.accession']))
    Mat1toMat2Overlap=length(overlap(Mat1Anchor[,'protein.accession'],Mat2))
    #draw.pairwise.venn(Mat1Area,Mat2Area,Mat1toMat2Overlap,category=c("ORFs","Uniprot"),euler.d=TRUE,scaled=TRUE,ext.text = TRUE,fill=c("black","white"))
    print(paste("Mat1 Anchor protein:",length(unique(Mat1Anchor[,'protein.accession']))))
    print(paste("Mat2 Anchor protein:",length(unique(Mat2Anchor[,'protein.accession']))))
    print(paste("Mat1 All protein:",length(unique(Mat1[,'protein.accession']))))
    print(paste("Mat2 All protein:",length(unique(Mat2[,'protein.accession']))))
    print(paste("Overlap: Mat1 to Mat2:",length(overlap(Mat1Anchor[,'protein.accession'],Mat2))))
    print(paste("Overlap: Mat2 to Mat1:",length(overlap(Mat2Anchor[,'protein.accession'],Mat1))))
    print(paste("Overlap:Anchor to Anchor",length(overlap(Mat1Anchor[,'protein.accession'],Mat2Anchor))))
    print(paste("Overlap:All to All",length(overlap(Mat1[,'protein.accession'],Mat2))))
    Mat1Anchor_no_iso = voidIsoformEffect(Mat1Anchor[,'protein.accession'])
    Mat2Anchor_no_iso = voidIsoformEffect(Mat2Anchor[,'protein.accession'])
    
    print(paste("Mat1 PAG ignoring isoform:",length(unique(Mat1Anchor_no_iso))))
    print(paste("Mat2 PAG ignoring isoform:",length(unique(Mat2Anchor_no_iso))))
    
    print(paste("Overlap without isoform:", length(overlap2(unique(Mat1Anchor_no_iso),unique(Mat2Anchor_no_iso)))))
    overlapUniPAG=unique(Mat2[overlapIndices2(Mat2[,'protein.accession'], unique(Mat1[,'protein.accession'])),'PAG.ID'])
    print(paste("Overlap:Uniprot all to Trinity anchor",length(overlap2(Mat2[,'protein.accession'], unique(Mat1Anchor[,'protein.accession'])))))
    print(paste("Overlap:Uniprot all to Trinity anchor PAGs",length(overlapUniPAG)))
    print(paste("Only:Uniprot all to Trinity anchor no of PAG",length(unique(Mat2[which(! Mat2[,'PAG.ID'] %in% overlapUniPAG),'PAG.ID']))))
    
    overlapTrPAG=unique(Mat1[overlapIndices2(Mat1[,'protein.accession'], unique(Mat2[,'protein.accession'])),'PAG.ID'])
    print(paste("Overlap:Trinity all to Uniprot anchor",length(overlap2(Mat1[,'protein.accession'], unique(Mat2Anchor[,'protein.accession'])))))
    print(paste("Overlap:Trinity all to Uniprot anchor PAGs",length(overlapTrPAG)))
    print(paste("Trinity only PAG group:",length(unique(Mat1[which(! Mat1[,'PAG.ID' %in% overlapTrPAG]),'PAG.ID']))))
    
    print(paste("Only:Trinity all to Uniprot anchor no of PAG",length(unique(Mat1[overlapIndices2(Mat1[,'protein.accession'], unique(Mat2[,'protein.accession'])),'PAG.ID']))))
    
    #write.table(Mat2[which(!Mat2[,'PAG.ID'] %in% Mat2[overlapIndices2(Mat2[,'protein.accession'], unique(Mat1[,'protein.accession'])),'PAG.ID']),],file=paste(d1,"only_uni_PAG",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat1[which(!Mat1[,'PAG.ID'] %in% Mat1[overlapIndices2(Mat1[,'protein.accession'], unique(Mat2[,'protein.accession'])),'PAG.ID']),],file=paste(d1,"only_ORFsT_PAG",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat1[overlapIndices(Mat1[,'protein.accession'],Mat2),],file=paste(d1,"overlap_all_",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat1Anchor[overlapIndices(Mat1Anchor[,'protein.accession'],Mat2),],file=paste(d1,"overlap",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat1Anchor[exclusive(Mat1Anchor[,'protein.accession'],Mat2),],file=paste(d1,"only",f1,upper,"_compared_",f2,"_",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat2Anchor[exclusive(Mat2Anchor[,'protein.accession'],Mat1),],file=paste(d2,"only",f2,upper,"_compared_",f1,"_",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #list(Mat1_back,Mat2_back)
}

########################################################################################################

##Species finding
length(grep("_HUMAN$",unique(TAdenoFAnchor[overlapIndices(TAdenoFAnchor[,'protein.accession'],UAdenoFAnchor),'protein.accession'])))
dir1="C:/Users/shyama/Dropbox/results/Human_adenovirus/"
dir2="C:/Users/shyama/Dropbox/results/Human_adenovirus/"
filename1="DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping+prt.csv"
filename2="DM_from_raw_trinity_uniprot+fdr+th+grouping+prt.csv"
filename3="blast_trinity_id_clean.csv"

Mat1=proteinGroup(filename1,dir1)
Mat2=proteinGroup(filename2,dir1)

##might need to manipulate 'protein.accession' values.

rev=1
peptide=1
pepThreshold=1
Mat1Filtered=proteinGroupFiltered(proteinGroup(filename1,dir1),rev,peptide,pepThreshold)
Mat2Filtered=proteinGroupFiltered(proteinGroup(filename2,dir1),rev,peptide,pepThreshold)
Mat1FAnchor=proteinGrpAnchor(Mat1Filtered)
Mat2FAnchor=proteinGrpAnchor(Mat2Filtered)
################Length###################
UT2_length=as.matrix(read.csv(file=paste(dir1,"trinity_uniprot_no_duplicatev2.tabular",sep=""),header=TRUE))
UT1_length=as.matrix(read.csv(file=paste(dir1,"trinity_uniprotv2.tabular",sep=""),header=TRUE))
Mat1Filtered[,'protein.accession']=sub("sw","sp",Mat1Filtered[,'protein.accession'])
Mat2Filtered[,'protein.accession']=sub("sw","sp",Mat2Filtered[,'protein.accession'])
lengthUT1Found=UT1_length[which(UT1_length[,1] %in% Mat2Filtered[,'protein.accession']),]
lengthUT2Found=UT2_length[which(UT2_length[,1] %in% Mat1Filtered[,'protein.accession']),]
lengthUT1notFound=UT1_length[which(!UT1_length[,1] %in% Mat2Filtered[,'protein.accession']),]
lengthUT2notFound=UT2_length[which(!UT2_length[,1] %in% Mat1Filtered[,'protein.accession']),]

 write.table(lengthUT2notFound,file=paste(folder,"humanAdenoUT2UnIdentifiedProteinLength.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
> write.table(lengthUT1notFound,file=paste(folder,"humanAdenoUT1UnIdentifiedProteinLength.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
###################################
blastMat=blast(filename3,dir2)
flag=1
blastMatFiltered=blastFilterMatch(blastMat,flag)
upper=0.1
lower=-0.9
blastMatFEval=blastFilterEval(blastMatFiltered,upper)
##might need to manipulate 'query_name' and 'hit_def' column
column=4
firstWord(blastMatFEval,column)
rIds=replaceIds(Mat1FAnchor,blastMatFEval)
overlap(rIds,Mat2FAnchor)
exclusive(rIds,Mat2FAnchor)

draw.pairwise.venn(2921, 2919, 2913, category = rep("", 2), euler.d = TRUE,
scaled = TRUE, inverted = FALSE, ext.text = TRUE, ext.percent = rep(0.05, 3),
lwd = rep(2, 2), lty = rep("solid", 2), col = rep("black", 2), fill = NULL,
alpha = rep(0.5, 2), label.col = rep("black", 3), cex = rep(1, 3),
fontface = rep("plain", 3), fontfamily = rep("serif", 3), cat.pos = c(-50, 50),
cat.dist = rep(0.025, 2), cat.cex = rep(1, 2), cat.col = rep("black", 2),
cat.fontface = rep("plain", 2), cat.fontfamily = rep("serif", 2),
cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos = "outer",
cat.prompts = FALSE,
ext.pos = rep(0, 2), ext.dist = rep(0, 2), ext.line.lty = "solid",
ext.length = rep(0.95, 2), ext.line.lwd = 1, rotation.degree = 0,
rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05, offset = 0)

#############################  identified protein length  #######################################

blastMatgetORF[which(blastMatgetORF[,'query_name'] %in% MatTAdenoF[,'protein.accession']),'query_length']

write.table(lengthAdenoTNot,file=paste(folder,"humanAdenogetORFUnidentifiedProteinLength.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
write.table(lengthAdenoTPITNot,file=paste(folder,"humanAdenoPITORFUnidentifiedProteinLength.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
write.table(lengthAdenoUNot,file=paste(folder,"humanAdenoUniUnidentifiedProteinLength.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)

rev=1
peptide=1
pepThreshold=1
upper=0.1
lower=-0.9

##########################Bat##################################
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



############################# length boxplot ########################################
 library(reshape)
 library(ggplot2)
a = read.delim("C:/Users/shyama/Dropbox/results/Human_adenovirus/humanAdenoProteinLengths.csv",sep=",",header=T)

X = rep( 1:dim(a)[1], 10)
a_melted = cbind( melt(a) , X )


## Line Plot ##
p1 <- ggplot(a_melted, aes( x = X, y = value, colour=variable, group=variable)) + geom_line() + ggtitle(" Database Length")
p1
## Box Plot ##
p2 <- ggplot(a_melted, aes( factor( variable ), value, colour=variable ) ) +xlab("Database") + ylab("Length") + geom_boxplot() + ggtitle("Protein/ORFs Length")
p2
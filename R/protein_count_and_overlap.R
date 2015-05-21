

########################################################################################################
###############                     Overlap of two MSGF+ results set                 ###################                            
########################################################################################################

### Parameter description ###
### f1 is proteinGroup file name ###
### f2 is proteinGroup file name ###
### f3 is blastCSV file name ###
### d1 is folder location of f1 ###
### d2 is folder location of f2 ###
### d3 is folder location of f3 ###
### Make sure f1 result database is the query for blast ###

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

##### create function for these#####

PAG_overlap<-function(M1,M2)
{
    count=0
    M1_temp=M1
    print(dim(M1))
    idxList=c()
    for(i in 1:dim(M1)[1])
    {
        indx=which(M2[,'protein.accession']==M1[i,'protein.accession'])
        
        if(length(indx)>0)
        {
            count=count+1
            idxList=c(idxList,i)
            M2=M2[-which(M2[,'PAG.ID']==M2[indx,'PAG.ID']),]
        }
    }
    M1_temp=M1[-idxList,]
    print("M1_temp")
    print(dim(M1_temp))
    list(count=count,M1=M1_temp,M2=M2,idx=idxList)
}
 m1=Mat1Anchor[!duplicated(Mat1Anchor[,'protein.accession']),]
 Mat1Sub=as.matrix(Mat1[grep('sequence.*',Mat1[,'group.membership']),]) #Mat1 sub members
 #removing any duplicate protein, if it is uniprot, then there should not be any duplicate. if it is PIT result, there will be duplicates because often multiple ORFs map back to same uniprot protein.
 m1_sub=Mat1Sub[!duplicated(Mat1Sub[,'protein.accession']),]
 #finding submember protein ids that does not exist in the anchor protein list.
m1_sub_dup=m1_sub[which(!m1_sub[,'protein.accession'] %in% m1[,'protein.accession']),]
#findinng PAGs that only has submembers. 
m1_sub_dup_no_anchor=m1_sub_dup[which(!m1_sub_dup[,'PAG.ID'] %in% m1[,'PAG.ID']),]
#removing submembers without an anchor protein.
m1_sub_dup_anchor=m1_sub_dup[-which(!m1_sub_dup[,'PAG.ID'] %in% m1[,'PAG.ID']),]
m1_merged=m1
#merged unique protein/PAg lists. duplicates from anchor and submembers are done separately to avoid removing an anchor id when there is a submember with same id. 
m1_merged=rbind(m1_merged, m1_sub_dup_anchor)


a_to_a=m1[overlap2(m1[,'protein.accession'],Mat2Anchor[,'protein.accession']),]
 
 Mat2atoaPAG=Mat2Anchor[overlap2(Mat2Anchor[,'protein.accession'],m1[,'protein.accession']),'PAG.ID']
 Mat2_without_atoaPAG=Mat2[-which(Mat2[,'PAG.ID'] %in% Mat2atoaPAG),]
 Mat2_without_atoaPAG_anchor=proteinGrpAnchor(Mat2_without_atoaPAG)

 
  m1_merge_atoaPAG=a_to_a[,'PAG.ID']
 m1_merged_without_atoaPAG=m1_merged[-which(m1_merged[,'PAG.ID'] %in% m1_merge_atoaPAG),]
 
 m1_without_atoaPAG_anchor=proteinGrpAnchor(m1_merged_without_atoaPAG)
 a_to_mapped_pags=m1_merged[which(m1_merged[,'PAG.ID'] %in% a_to_a[,'PAG.ID']),]
 
 ##anchor to anchor maps and their protein group members are written to a file. But overlap count is only those are anchor proteins.
 write.table(a_to_mapped_pags,file=paste(d1,"AnchortoAnchor_PAGs_uni_pasa",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
 
 res=PAG_overlap(m1_without_atoaPAG_anchor,Mat2_without_atoaPAG)
 trinityAnchortoUniSubOverlap=m1_merged_without_atoaPAG[which(m1_merged_without_atoaPAG[,'PAG.ID'] %in% m1_without_atoaPAG_anchor[res$idx,'PAG.ID']),]
 #write.table(trinityAnchortoUniSubOverlap,file=paste(d1,"PASAAnchortoUniSub_PAGs",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
  dim(a_to_a)
 print(res$count)
 m1_merged_without_atosub=m1_merged_without_atoaPAG[which(m1_merged_without_atoaPAG[,'PAG.ID'] %in% res$M1[,'PAG.ID']),]
 dim(m1_merged_without_atosub)
 dim(res$M1)
 dim(res$M2)
 Mat2_without_subtoa_anchor = proteinGrpAnchor(res$M2)
 dim(Mat2_without_subtoa_anchor)
 
 #m1_without_atoaPAG_anchor[res$idx,'protein.accession']
 res2=PAG_overlap(Mat2_without_subtoa_anchor,m1_merged_without_atosub)
 print(res2$count)
 dim(res2$M1)
 dim(proteinGrpAnchor(res2$M2))
 res2$idx
 Mat2_without_subtoa_anchor[res2$idx,]
 idx=which(m1_merged_without_atosub[,'protein.accession'] %in% Mat2_without_subtoa_anchor[res2$idx,'protein.accession'])
  uniAnchortoTrinityOverlap=m1_merged_without_atosub[which(m1_merged_without_atosub[,'PAG.ID'] %in% m1_merged_without_atosub[idx,'PAG.ID']),]
#write.table( uniAnchortoTrinityOverlap ,file=paste(d1,"UniAnchortoCuffSub_PAGs_uni_trinity",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
 anchorToSubOverlap=rbind(trinityAnchortoUniSubOverlap,uniAnchortoTrinityOverlap)
 write.table( anchorToSubOverlap ,file=paste(d1,"AnchortoSub_PAGs_uni_PASA",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
 
 overlapMat=rbind(a_to_mapped_pags,anchorToSubOverlap)
 write.table( overlapMat ,file=paste(d1,"overlap_PAGs_uni_PASA",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
 write.table(res2$M1,file=paste(d1,"PAGComparison_UniprotOnly_withPASA",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
 write.table(proteinGrpAnchor(res2$M2),file=paste(d1,"PAGComparison_PASAOnly_withUniprot",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)

cuffOnly=proteinGrpAnchor(res2$M2)
cuffOnlyIds=cuffOnly[grep("^CUFF",cuffOnly[,'protein.accession']),'protein.accession']
cuffNovelMap=blastDB[which(blastDB[,'query_name'] %in% cuffOnlyIds),]


 #length(overlap2(m1_without_atoaPAG_anchor[,'protein.accession'], Mat2_without_atoaPAG[,'protein.accession']))
 #m1_atosubPAG=m1_without_atoaPAG_anchor[which(m1_without_atoaPAG_anchor[,'protein.accession'] %in% Mat2_without_atoaPAG[,'protein.accession']),'PAG.ID']
 
 #m1_merged_without_atosubPAG=m1_merged_without_atoaPAG[-which(m1_merged_without_atoaPAG[,'PAG.ID'] %in% unique(m1_atosubPAG)),]
#dim(proteinGrpAnchor(m1_merged_without_atosubPAG))
 
 #Mat2_subtoAPAGs=Mat2_without_atoaPAG[which(Mat2_without_atoaPAG[,'protein.accession'] %in% unique(m1_without_atoaPAG_anchor[,'protein.accession'])),'PAG.ID']
 #Mat2_without_atosubPAG=Mat2_without_atoaPAG_anchor[-which(Mat2_without_atoaPAG_anchor[,'PAG.ID'] %in% unique(Mat2_subtoAPAGs)),]
 #Mat2_without_atosubPAGAnchor=proteinGrpAnchor(Mat2_without_atosubPAG)

 #length(unique(Mat2_subtoAPAGs))
 #length(overlap2(Mat2_without_atosubPAGAnchor[,'protein.accession'],m1_merged_without_atosubPAG[,'protein.accession']))
########## Venn Diagram ########
library(VennDiagram)
draw.pairwise.venn(Mat1Area,Mat2Area,Mat1toMat2Overlap,category=c("ORFs","Uniprot"),euler.d=TRUE,scaled=TRUE,ext.text = TRUE,fill=c("black","white"))

venn.plot <- venn.diagram(list(Uniprot = 1:26401, ORFs = 1884:27789),"Venn_2set_PeptideTrinity.tiff",fill=c("black","white"), alpha=c(0.65,.75), offset=0.5,euler.d=TRUE,scaled=TRUE, ext.text=TRUE, cex=c(2,2,2), cat.col=c("white","white"), ext.line.lwd=0, label.col=c("white","black","white"), ext.percent=c(.25,.25,.25))
venn.plot <- venn.diagram(list(Uniprot = 1:25970, ORFs = 5216:26816),"Venn_2set_PeptideCuff.tiff",fill=c("black","white"), alpha=c(0.65,.75), offset=0.5,euler.d=TRUE,scaled=TRUE, ext.text=TRUE, cex=c(2,2,2), cat.col=c("white","white"), ext.line.lwd=0, label.col=c("white","black","white"), ext.percent=c(.25,.25,.25))
venn.plot <- venn.diagram(list(Uniprot = 1:2986, ORFs = 332:3041),"Venn_2set_ProteinCuff.tiff",fill=c("black","white"), alpha=c(0.65,.75), offset=0.5,euler.d=TRUE,scaled=TRUE, ext.text=TRUE, cex=c(2,2,2), cat.col=c("white","white"), ext.line.lwd=0, label.col=c("white","black","white"), ext.percent=c(.25,.25,.25))

venn.plot <- venn.diagram(list(Uniprot = 1:3015, ORFs = 217:3143),"Venn_2set_Protein.tiff",fill=c("black","white"), alpha=c(0.65,.75), offset=0.5,euler.d=TRUE,scaled=TRUE, ext.text=TRUE, cex=c(2,2,2), cat.col=c("white","white"), ext.line.lwd=0, label.col=c("white","black","white"), ext.percent=c(.25,.25,.25))

venn.plot <- venn.diagram(
x = list(
A = c(1:90),
B = c(11:90),
C = c(81:90)
),
euler.d = TRUE,
filename = "Euler_3set_scaled.tiff",
cex = 2.5,
cat.cex = 2.5,
cat.pos = 0
);

########################################## 18-05-2015, PASA uniprot overlap #########################################################
f1="pasa_assemblyV1+fdr+th+grouping+prt.csv"
f2="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
f3="human_adeno_mydb_pasa.assemblies_ORFsV1.csv"
d1="D:/data/Results/Human-Adeno/Identification/PASA/sORF/"
d2="D:/data/Results/Human-Adeno/Identification/"
#d1=d2
d3="D:/data/blast/blastCSV/"
upper=0.000000000000000000000000000001


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


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
### Make sure f1 protein accessions are the query for blast ###

source("D:/Code/Proteomics/R/RLib.R")

readMatsAndPostProcess <- function(f1,f2,f3,d1,d2,d3,rev=1,peptide=1,pepThreshold=1,upper=0.1)
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
    ##Unique proteins after replacing Mat1 ids with Mat2
    print(paste("Mat1 unique protein count after replacement:",length(unique(Mat1[,'protein.accession']))))
    list(Mat1=Mat1, Mat2=Mat2)
}

readMatsAndPostProcessBothMappedToThirdDB <- function(f1,f2,f3,f4,d1,d2,d3,d4,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    Mat2=proteinGroupFiltered(proteinGroup(f2,d2),rev,peptide,pepThreshold)
    blastDB1=blastFilterEval(blastFilterMatch(blast(f3,d3),1),upper)
    blastDB2=blastFilterEval(blastFilterMatch(blast(f4,d4),1),upper)
    Mat1=replaceComma(Mat1,3)
    Mat1=upto(Mat1,3,';')
    Mat1=upto(Mat1,3,' ')
    Mat2=replaceComma(Mat2,3)
    Mat2=upto(Mat2,3,';')
    Mat2=upto(Mat2,3,' ')
    blastDB1=upto(blastDB1,1,';')
    blastDB1=upto(blastDB1,1,',')
    blastDB1=upto(blastDB1,4,' ')
    
    blastDB2=upto(blastDB2,1,';')
    blastDB2=upto(blastDB2,1,',')
    blastDB2=upto(blastDB2,4,' ')
    
    Mat1[,'protein.accession']=sub("^sw","sp",Mat1[,'protein.accession'])
    Mat2[,'protein.accession']=sub("^sw","sp",Mat2[,'protein.accession'])
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB1)
    Mat2[,'protein.accession']=replaceIds(Mat2,blastDB2)
    

    ##Unique proteins after replacing Mat1 ids with Third DB
    print(paste("Mat1 unique protein count after replacement:",length(unique(Mat1[,'protein.accession']))))
    ##Unique proteins after replacing Mat2 ids with Third DB
    print(paste("Mat2 unique protein count after replacement:",length(unique(Mat2[,'protein.accession']))))
    
    list(Mat1=Mat1,Mat2=Mat2)
}

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

calculateAndWrite <- function(Mat1, Mat2, Mat1FileNamePostfix, Mat2FileNamePostfix)
{
    Mat1Anchor=proteinGrpAnchor(Mat1)
    Mat2Anchor=proteinGrpAnchor(Mat2)
    
    ## After reaplacing Mat1 with Mat2 ids, Mat1 ends up with non-unique protein accessions. In next line is removing duplicates from anchor protein ids.
    m1=Mat1Anchor[!duplicated(Mat1Anchor[,'protein.accession']),]
    print(paste("m1:",dim(m1)))
    Mat1Sub=as.matrix(Mat1[grep('sequence.*',Mat1[,'group.membership']),]) #Mat1 sub members
    #removing any duplicate protein, if it is uniprot, then there should not be any duplicate. if it is PIT result, there will be duplicates because often multiple ORFs map back to same uniprot protein.
    m1_sub=Mat1Sub[!duplicated(Mat1Sub[,'protein.accession']),]
    #finding submember protein ids that does not exist in the anchor protein list.
    m1_sub_dup=m1_sub[which(!m1_sub[,'protein.accession'] %in% m1[,'protein.accession']),]
    #findinng PAGs that only has submembers.
    m1_sub_dup_no_anchor=m1_sub_dup[which(!m1_sub_dup[,'PAG.ID'] %in% m1[,'PAG.ID']),]
    #removing submembers without an anchor protein.
    m1_sub_dup_anchor=m1_sub_dup[-which(!m1_sub_dup[,'PAG.ID'] %in% m1[,'PAG.ID']),]
    print(paste("m1_sub_dup_anchor:",dim(m1_sub_dup_anchor)))
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
    write.table(a_to_mapped_pags,file=paste(d1,"AnchortoAnchor_PAGs_",Mat1FileNamePostfix,"_",Mat2FileNamePostfix,upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    
    res=PAG_overlap(m1_without_atoaPAG_anchor,Mat2_without_atoaPAG)
    trinityAnchortoUniSubOverlap=m1_merged_without_atoaPAG[which(m1_merged_without_atoaPAG[,'PAG.ID'] %in% m1_without_atoaPAG_anchor[res$idx,'PAG.ID']),]
    #write.table(trinityAnchortoUniSubOverlap,file=paste(d1,"PASAAnchortoUniSub_PAGs",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    
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
    anchorToSubOverlap=rbind(trinityAnchortoUniSubOverlap,uniAnchortoTrinityOverlap)
    write.table( anchorToSubOverlap ,file=paste(d1,"AnchortoSub_PAGs_",Mat1FileNamePostfix,"_",Mat2FileNamePostfix,upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
    
    overlapMat=rbind(a_to_mapped_pags,anchorToSubOverlap)
    write.table( overlapMat ,file=paste(d1,"overlap_PAGs_",Mat1FileNamePostfix,"_",Mat2FileNamePostfix,upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
    write.table(res2$M1,file=paste(d1,"PAGComparison_",Mat2FileNamePostfix,"Only_with",Mat1FileNamePostfix,upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
    write.table(proteinGrpAnchor(res2$M2),file=paste(d1,"PAGComparison_",Mat1FileNamePostfix,"Only_with",Mat2FileNamePostfix,upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
    print("a_to_a")
    print(dim(a_to_a))
    print("res$count")
    print(res$count)
    print("res2$count")
    print(res2$count)
}

calculateAndReturnOverlapMatrix <- function(Mat1, Mat2, Mat1FileNamePostfix, Mat2FileNamePostfix)
{
    Mat1Anchor=proteinGrpAnchor(Mat1)
    Mat2Anchor=proteinGrpAnchor(Mat2)
    
    ## After reaplacing Mat1 with Mat2 ids, Mat1 ends up with non-unique protein accessions. In next line is removing duplicates from anchor protein ids.
    m1=Mat1Anchor[!duplicated(Mat1Anchor[,'protein.accession']),]
    print(paste("m1:",dim(m1)))
    Mat1Sub=as.matrix(Mat1[grep('sequence.*',Mat1[,'group.membership']),]) #Mat1 sub members
    #removing any duplicate protein, if it is uniprot, then there should not be any duplicate. if it is PIT result, there will be duplicates because often multiple ORFs map back to same uniprot protein.
    m1_sub=Mat1Sub[!duplicated(Mat1Sub[,'protein.accession']),]
    #finding submember protein ids that does not exist in the anchor protein list.
    m1_sub_dup=m1_sub[which(!m1_sub[,'protein.accession'] %in% m1[,'protein.accession']),]
    #findinng PAGs that only has submembers.
    m1_sub_dup_no_anchor=m1_sub_dup[which(!m1_sub_dup[,'PAG.ID'] %in% m1[,'PAG.ID']),]
    #removing submembers without an anchor protein.
    m1_sub_dup_anchor=m1_sub_dup[-which(!m1_sub_dup[,'PAG.ID'] %in% m1[,'PAG.ID']),]
    print(paste("m1_sub_dup_anchor:",dim(m1_sub_dup_anchor)))
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
    
    
    res=PAG_overlap(m1_without_atoaPAG_anchor,Mat2_without_atoaPAG)
    trinityAnchortoUniSubOverlap=m1_merged_without_atoaPAG[which(m1_merged_without_atoaPAG[,'PAG.ID'] %in% m1_without_atoaPAG_anchor[res$idx,'PAG.ID']),]
    #write.table(trinityAnchortoUniSubOverlap,file=paste(d1,"PASAAnchortoUniSub_PAGs",upper,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    
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
    anchorToSubOverlap=rbind(trinityAnchortoUniSubOverlap,uniAnchortoTrinityOverlap)
    
    overlapMat=rbind(a_to_mapped_pags,anchorToSubOverlap)
    print("a_to_a")
    print(dim(a_to_a))
    print("res$count")
    print(res$count)
    print("res2$count")
    print(res2$count)
}

PAG_overlap_Main <- function(f1,f2,f3,d1,d2,d3,rev=1,peptide=1,pepThreshold=1,upper=0.1, Mat1FileNamePostfix, Mat2FileNamePostfix)
{
    Mat=readMatsAndPostProcess(f1,f2,f3,d1,d2,d3,rev,peptide,pepThreshold,upper)
    calculateAndWrite(Mat$Mat1, Mat$Mat2, Mat1FileNamePostfix, Mat2FileNamePostfix)
}

## use case: when results of X DB search and Y DB search is being compared but there is not mapping available between X and Y. Though Mapping
## between X and Z, and Y and Z is available.
PAG_overlap_ThirdDB_Map_Main <- function(f1,f2,f3,f4,d1,d2,d3,d4,rev=1,peptide=1,pepThreshold=1,upper=0.1, Mat1FileNamePostfix, Mat2FileNamePostfix)
{
    Mat=readMatsAndPostProcessBothMappedToThirdDB(f1,f2,f3,f4,d1,d2,d3,d4,rev,peptide,pepThreshold,upper)
    calculateAndWrite(Mat$Mat1, Mat$Mat2, Mat1FileNamePostfix, Mat2FileNamePostfix)
}



##################################################################################################################################################
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


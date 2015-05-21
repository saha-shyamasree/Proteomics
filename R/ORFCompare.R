### this code compare ORFs to Uniprot proteins blast result, where ORFs are generated from either trinity assebled transcripts or PASA assemblies.

> dir="D:/data/blast/blastCSV/"
> tFile="trinityV5.csv"
> pFile="human_adeno_mydb_pasa.assemblies_ORFs.csv"
upper=0.000000000000000000000000000001
trinityMat=blastFilterEval(blastFilterMatch(blast(tFile,dir),1),upper)
pasaMat=blastFilterEval(blastFilterMatch(blast(pFile,dir),1),upper)

tMat=blast(tFile,dir)
pMat=blast(pFile,dir)

library(reshape)
library(ggplot2)
###################################################################################################################################
LENGTH
###################################################################################################################################
#######uniprot homologous ############################################
pL=c(pasaMat[,'query_length'],rep(NA,dim(trinityMat)[1]-dim(pasaMat)[1]))
mergedL=cbind(pasa=pL,trinity=trinityMat[,'query_length'])
mergedL_melted = melt(mergedL)
colnames(mergedL_melted)=c('SerialNo','Software','Length')
## Line Plot ##
p1 <- ggplot(mergedL_melted, aes( x = Software, y = Length, colour=Software, group=Software)) + geom_line() + ggtitle(" Database Length")
p1
##box plot
p2 <- ggplot(mergedL_melted, aes( factor( Software ), Length, colour=Software ) ) +xlab("Database") + ylab("Length") + geom_boxplot() + ggtitle("ORFs Length: Trinity vs PASA (homologous to uniprot)")

##density
ggplot(mergedL_melted, aes(value, fill = Software)) + geom_density(alpha = 0.2)

#### all #############################################

pL=c(pMat[,'query_length'],rep(NA,dim(tMat)[1]-dim(pMat)[1]))
mergedL=cbind(pasa=pL,trinity=tMat[,'query_length'])
mergedL_melted = melt(mergedL)
colnames(mergedL_melted)=c('SerialNo','Software','Length')
## Line Plot ##
p1 <- ggplot(mergedL_melted, aes( x = Software, y = Length, colour=Software, group=Software)) + geom_line() + ggtitle(" Database Length")
p1
##box plot
p2 <- ggplot(mergedL_melted, aes( factor( Software ), Length, colour=Software ) ) +xlab("Database") + ylab("Length") + geom_boxplot() + ggtitle("ORFs Length: Trinity vs PASA (All)")

##density
ggplot(mergedL_melted, aes(value, fill = Software)) + geom_density(alpha = 0.2)

###################################################################################################################################
identities/hit_length
###################################################################################################################################
#######uniprot homologous ############################################
pL=c(pasaMat[,'long_match'],rep(NA,dim(trinityMat)[1]-dim(pasaMat)[1]))
mergedL=cbind(pasa=pL,trinity=trinityMat[,'long_match'])
mergedL_melted = melt(mergedL)
colnames(mergedL_melted)=c('SerialNo','Software','long_match')
##box plot
p2 <- ggplot(mergedL_melted, aes( factor( Software ), long_match, colour=Software ) ) +xlab("Software") + ylab("identities to hit length ration") + geom_boxplot() + ggtitle("ORFs identities to hit length ratio: Trinity vs PASA (homologous to uniprot)")

##density
ggplot(mergedL_melted, aes(long_match, fill = Software)) + geom_density(alpha = 0.2)


###################################################################################################################################
identities/alignment_length
###################################################################################################################################
#######uniprot homologous ############################################
pL=c(pasaMat[,'good_match'],rep(NA,dim(trinityMat)[1]-dim(pasaMat)[1]))
mergedL=cbind(pasa=pL,trinity=trinityMat[,'good_match'])
mergedL_melted = melt(mergedL)
colnames(mergedL_melted)=c('SerialNo','Software','good_match')
##box plot
p2 <- ggplot(mergedL_melted, aes( factor( Software ), good_match, colour=Software ) ) +xlab("Software") + ylab("identities to alignment length ratio") + geom_boxplot() + ggtitle("ORFs identities to alignment length ratio: Trinity vs PASA (homologous to uniprot)")

##density
ggplot(mergedL_melted, aes(good_match, fill = Software)) + geom_density(alpha = 0.2)


############# above threshold ##############################
th=0.9
tempP=pasaMat[which(pasaMat[,'good_match']>th),'good_match']
tempT=trinityMat[which(trinityMat[,'good_match']>th),'good_match']
pL=c(tempP,rep(NA,length(tempT)[1]-length(tempP)))
mergedL=cbind(pasa=pL,trinity=tempT)
mergedL_melted = melt(mergedL)
colnames(mergedL_melted)=c('SerialNo','Software','good_match')
##box plot
p2 <- ggplot(mergedL_melted, aes( factor( Software ), good_match, colour=Software ) ) +xlab("Software") + ylab("identities to alignment length ratio") + geom_boxplot() + ggtitle("ORFs identities to alignment length ratio: Trinity vs PASA (homologous to uniprot)")

##density
ggplot(mergedL_melted, aes(good_match, fill = Software)) + geom_density(alpha = 0.2)

##########################################################################################
Comparing protein identification results for ORFs from trinity assembly and PASA assembly.
##########################################################################################

pPrtFile="pasa_assemblyV1+fdr+th+grouping+prt.csv"
tPrtFile="trinity_PITORF+fdr+th+grouping+prtV5.csv"
pPrtDir="D:/data/Results/Human-Adeno/Identification/"
tPrtDir="C:/Users/shyama/Dropbox/results/Human_adenovirus/sORF/"
pBlastFile="human_adeno_mydb_pasa.assemblies_ORFs.csv"
tBlastFile="trinityV5.csv"
blastDir="D:/data/blast/blastCSV/"

upper=0.000000000000000000000000000001

#tBMat=blastFilterEval(blastFilterMatch(blast(tBlastFile,blastDir),1),upper)
#pBMat=blastFilterEval(blastFilterMatch(blast(pBlastFile,blastDir),1),upper)

matList=twoDataSetsOverlap(pPrtFile, tPrtFile, pBlastFile, tBlastFile, pPrtDir, tPrtDir, blastDir, blastDir, 1,1,1,upper)

> pAncSameSet=proteinGrpAnchorAndSameSet(matList[[1]])
> tAncSameSet=proteinGrpAnchorAndSameSet(matList[[2]])

##Mat1 holds anchor proteins and same set proteins
differentORFtosameUniFromSamePAG<-function(Mat1)
{
    anchors=proteinGrpAnchor(Mat1)
    sameSet=proteinGrpSameSet(Mat1)
    apply(anchors,1,function(x,y){ which(y[which(y[,'PAG.ID']==x['PAG.ID']),'protein.accession']==x['protein.accession'])}, y=sameSet)
}

twoDataSetsOverlap <- function(f1,f2,f3,f4,d1,d2,d3,d4,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    Mat2=proteinGroupFiltered(proteinGroup(f2,d2),rev,peptide,pepThreshold)
    blastDB1=blastFilterEval(blastFilterMatch(blast(f3,d3),1),upper)
    blastDB2=blastFilterEval(blastFilterMatch(blast(f4,d4),1),upper)
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
    blastDB1=upto(blastDB1,1,';')
    blastDB1=upto(blastDB1,1,',')
    blastDB1=upto(blastDB1,4,' ')
    
    blastDB2=upto(blastDB2,1,';')
    blastDB2=upto(blastDB2,1,',')
    blastDB2=upto(blastDB2,4,' ')
    Mat1[,'protein.accession']=sub("^sw","sp",Mat1[,'protein.accession'])
    print(dim(Mat1))
    Mat2[,'protein.accession']=sub("^sw","sp",Mat2[,'protein.accession'])
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB1)
    Mat2[,'protein.accession']=replaceIds(Mat2,blastDB2)
    list(Mat1,Mat2)
    #print(dim(Mat1))
    #print("Mat1 overlap dim")
    #print(dim(Mat1))
    #write.table(Mat1,file=paste(d1,f1,"tempMat1.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat2,file=paste(d2,f2,"tempMat2.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #print(length(overlapIndices(proteinGrpAnchor(Mat1)[,'protein.accession'],Mat2)))
    #Mat1Anchor=proteinGrpAnchor(Mat1)
    #Mat2Anchor=proteinGrpAnchor(Mat2)
}

overlap <- function()
{
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
}

blast<-function(filename,folder)
{
    read.csv(file=paste(folder,filename,sep=""),header=TRUE,quote = "")
}

blastFilterMatch <- function(blastMat,flag)
{
    if(flag==1)
    {
        blastMat[which(blastMat[,'match']=='yes'),]
    }
    else
    {
        if(flag==0)
        {
            blastMat[which(blastMat[,'match']=='no'),]
        }
        else
        {
            NULL;
        }
    }
}

blastFilterEval <- function(blastMat,upper)
{
    print(paste("evalue=",upper))
    subset(blastMat,e.value<=upper)
}


proteinGroup <- function(filename,dirc)
{
    as.matrix(read.csv(file=paste(dirc,filename,sep=""),header=TRUE))
}

proteinGroupFiltered<- function(proteinGrp,rev,peptide,pepThreshold)
{
    proteinGrp=as.matrix(proteinGrp)
    print(dim(proteinGrp))
    if(rev==1)
    {
        proteinGrp=proteinGrp[-(grep("_REVERSED",proteinGrp[,'protein.accession'])),]
        print(dim(proteinGrp))
    }
    if(peptide==1)
    {
        proteinGrp[,'distinct.peptide.sequences']=as.numeric(proteinGrp[,'distinct.peptide.sequences'])
        proteinGrp=proteinGrp[which(proteinGrp[,'distinct.peptide.sequences']>pepThreshold),]
        print(dim(proteinGrp))
    }
    proteinGrp
}

proteinGrpSize <- function(proteinGrp)
{
    dim(unique(proteinGrp[,'PAG.ID']))
}

proteinGrpAnchor <- function(proteinGrp)
{
    as.matrix(proteinGrp[grep('anchor protein',proteinGrp[,'group.membership']),])
}

blast<-function(filename,folder)
{
    read.csv(file=paste(folder,filename,sep=""),header=TRUE)
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

blastFilterEval <- function(blastMat,upper,lower)
{
    #blastMat[,'e.value']=as.numeric(blastMat[,'e.value'])
    if(upper>=lower)
    {
        subset(blastMat,e.value<=upper)
    }
    else
    {
        NULL
    }
}


##separating id upto first space
firstWord<-function(Mat,column)
{
    Mat[,column]=sub("([^\\s]+) (.*)?","\\1",Mat[,column])
    Mat
}

upto<-function(Mat,column,character)
{
    Mat[,column]=sub(paste("([^",character,"]+)",character,"(.*)?",sep=""),"\\1",Mat[,column])
    Mat
}

##replacing any comma in a column value.
replaceComma<-function(Mat,column)
{
    Mat[,column]=gsub(",",";",Mat[,column])
    Mat
}

replaceIds<-function(Mat,B)
{
    ##B should hold blast result in CSV format
    ##Mat is protein identification result from MSGF+
    indices=which(Mat[,'protein.accession'] %in% B[,'query_name'])
    print("blast map count")
    print(length(indices))
    for(i in indices)
    {
        Mat[i,'protein.accession'] = as.character(B[which(B[,'query_name'] %in% Mat[i,'protein.accession']),'hit_def'])
    }
    Mat[,'protein.accession']
}

overlap<-function(replacedIds,Mat)
{
    which(unique(replacedIds) %in% Mat[,'protein.accession'])
}

overlapIndices<-function(replacedIds,Mat)
{
    which(replacedIds %in% Mat[,'protein.accession'])
}


#####non-anchor match analysis, mat1 should be the subset of original Mat1, where member of mat1 are anchor proteins from Mat1 that matches non-anchor proteins of Mat2
peptideCountCheck<-function(Mat1,Mat2)
{
    pag_ids=which(Mat2[,'description'] %in% Mat1[,'description'])
    indices=c()
    count=0
    for(i in 1:length(pag_ids))
    {
        pep=Mat2[Mat2[,'PAG.ID']==Mat2[pag_ids[i],'PAG.ID'] & Mat2[,'group.membership']=='anchor protein',c('protein.accession','distinct.peptide.sequences')]
        if(Mat2[pag_ids[i],'distinct.peptide.sequences']==pep[2])
        {
            indices=c(indices,pag_ids[i])
            print(paste("anchor pep count:",pep[2],",id",pep[1],",match pep count:",Mat2[pag_ids[i],'distinct.peptide.sequences'],",id:",Mat2[pag_ids[i],'protein.accession']))
            if(length(grep("-",pep[1]))>0 || length(grep("-",Mat2[pag_ids[i],'protein.accession']))>0)
            {
                count=count+1
            }
        }
    }
    print(count)
    Mat2[indices,'distinct.peptide.sequences']
}

exclusive<-function(replacedIds,Mat)
{
    which(!replacedIds %in% Mat[,'protein.accession'])
}

#B is blast result, Mat is protein group matrix from MSGF+, where Mat is similar to blast's column matches the Mat (i.e. if blastdb was creadted for uniprot proteins and trinity ORFs were blasted against it, then Mat should be anchor proteins.accession of protein groups exclusive to trinity ORFs database search)
exclusiveProteinLength<-function(Mat,B,column)
{
    which(unique(B[,column]) %in% Mat)
}
normalizedPAGScore<-function(Mat)
{
    pag=c()
    for(i in 1:dim(Mat)[1])
    {
        pag=c(pag,(as.numeric(Mat[i,'PAG.score'])/as.numeric(Mat[i,'distinct.peptide.sequences'])))
    }
    pag
}

proteinCountPerGroup<-function(Mat)
{
    anchor=proteinGrpAnchor(Mat)
    pag_prt_count=matrix(nrow=dim(anchor)[1],ncol=2)
    colnames(pag_prt_count)=c('PAG.ID','Protein.Count')
    for(i in 1:dim(anchor)[1])
    {
        pag_prt_count[i,'PAG.ID'] = anchor[i,'PAG.ID']
        pag_prt_count[i,'Protein.Count'] = length(which(Mat[,'PAG.ID']==anchor[i,'PAG.ID']))
    }
    pag_prt_count
}

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

twoDataSetsOverlap <- function(f1,f2,f3,d1,d2,d3,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    Mat1=proteinGroup(f1,d1)#proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    Mat2=proteinGroup(f2,d1)#proteinGroupFiltered(proteinGroup(f2,d1),rev,peptide,pepThreshold)
    blastDB=blastFilterEval(blastFilterMatch(blast(f3,d3),1),upper,0)
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
    Mat1_back=Mat1
    Mat2_back=Mat2
    Mat1=proteinGroupFiltered(Mat1_back,rev,peptide,pepThreshold)
    Mat2=proteinGroupFiltered(Mat2_back,rev,peptide,pepThreshold)
    #print(dim(Mat1))
    #print("Mat1 overlap dim")
    #print(dim(Mat1))
    #write.table(Mat1,file=paste(d1,"tempMat1.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #write.table(Mat2,file=paste(d2,"tempMat2.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    #print(length(overlapIndices(proteinGrpAnchor(Mat1)[,'protein.accession'],Mat2)))
    Mat1Anchor=proteinGrpAnchor(Mat1)
    Mat2Anchor=proteinGrpAnchor(Mat2)
    write.table(Mat1Anchor[overlapIndices(Mat1Anchor[,'protein.accession'],Mat2),],file=paste(d1,"overlap.tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    write.table(Mat1Anchor[exclusive(Mat1Anchor[,'protein.accession'],Mat2),],file=paste(d1,"only",f1,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    write.table(Mat2Anchor[exclusive(Mat2Anchor[,'protein.accession'],Mat1),],file=paste(d2,"only",f2,".tsv",sep=""),sep='\t',quote = FALSE,row.names = FALSE)
    list(Mat1_back,Mat2_back)
}

filenameBHT="trinity_bat+fdr+th+grouping+prt.csv"
filenameBHU="uniprot_bat+fdr+th+grouping+prt.csv"
folder="C:/Users/shyama/Dropbox/results/Bat_human_hendra/Bat/"
TBH=proteinGroupFiltered(proteinGroup(filenameBHT,folder),rev,peptide,pepThreshold)
UBH=proteinGroupFiltered(proteinGroup(filenameBHU,folder),rev,peptide,pepThreshold)
filenameDB="uni_db_p_alecto_blast.xml2.csv"
dbfolder="C:/Users/shyama/Dropbox/results/blastdb/blastCSV/"
blastBH=blastFilterEval(blastFilterMatch(blast(filenameDB,dbfolder),1),upper,lower)
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
blastMatFEval=blastFilterEval(blastMatFiltered,upper,lower)
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
blastBH=blastFilterEval(blastFilterMatch(blast(filenameDB,dbfolder),1),upper,lower)
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
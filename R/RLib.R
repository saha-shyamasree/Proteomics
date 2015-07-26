###PSMs and Peptides
readPSM<-function(filename,dir)
{
    as.matrix(read.csv(file=paste(dir,filename,sep=""),header=TRUE))
}

countPSMs<-function(Mat1)
{
    length(unique(Mat1[,'Spectrum.ID']))
}

overlapPepIndices<-function(Mat1, Mat2, column1, column2)
{
    which(Mat1[,column1] %in% Mat2[,column2])
}
overlapPep<-function(Mat1,Mat2,column1,column2)
{
    which(unique(Mat1[,column1]) %in% unique(Mat2[,column2]))
}
overlapSpectra<-function(Mat1,Mat2)
{
    length(which(unique(Mat1[,'Spectrum.ID']) %in% unique(Mat2[,'Spectrum.ID'])))
}

thresholdPassedPSMs<-function(Mat1)
{
    Mat1[which(Mat1[,'Pass.Threshold']=="TRUE"),]
}

belowThresholdPSMs<-function(Mat1)
{
    Mat1[which(Mat1[,'Pass.Threshold']=="FALSE"),]
}

####Proteins

proteinGroup <- function(filename,dirc)
{
    as.matrix(read.csv(file=paste(dirc,filename,sep=""),header=TRUE))
}
proteinList <- function(filename,dirc, sep)
{
    as.matrix(read.table(file=paste(dirc,filename,sep=""),sep=sep,header=TRUE))
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

oneHitWonderProteinGrp<-function(proteinGrp,rev)
{
    proteinGrp=as.matrix(proteinGrp)
    print(dim(proteinGrp))
    if(rev==1)
    {
        proteinGrp=proteinGrp[-(grep("_REVERSED",proteinGrp[,'protein.accession'])),]
        print(dim(proteinGrp))
    }
    proteinGrp[,'distinct.peptide.sequences']=as.numeric(proteinGrp[,'distinct.peptide.sequences'])
    proteinGrp=proteinGrp[which(proteinGrp[,'distinct.peptide.sequences']==1),]
}

proteinGrpSize <- function(proteinGrp)
{
    dim(unique(proteinGrp[,'PAG.ID']))
}

proteinGrpAnchor <- function(proteinGrp)
{
    as.matrix(proteinGrp[grep('anchor protein',proteinGrp[,'group.membership']),])
}

proteinGrpAnchorAndSameSet <- function(proteinGrp)
{
    as.matrix(proteinGrp[grep('anchor protein|sequence same-set protein',proteinGrp[,'group.membership']),])
}

proteinGrpSameSet <-function(proteinGrp)
{
    as.matrix(proteinGrp[grep('sequence same-set protein',proteinGrp[,'group.membership']),])
}

proteinGrpSub <- function(proteinGrp)
{
    as.matrix(proteinGrp[grep('sequence.*',proteinGrp[,'group.membership']),])
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

blastFilterEval <- function(blastMat,upper)
{
    print(paste("evalue=",upper))
    subset(blastMat,e.value<=upper)
}

blastFilterGoodMatch <- function(blastMat,threshold)
{
    #blastMat[,'e.value']=as.numeric(blastMat[,'e.value'])
    subset(blastMat,good_match>=threshold)
}

blastFilterLongMatch <- function(blastMat,threshold)
{
    #blastMat[,'e.value']=as.numeric(blastMat[,'e.value'])
    subset(blastMat,long_match>=threshold)
}

blastFilterLengthRatioMatch <- function(blastMat,threshold)
{
    if(!"length_ratio" %in% colnames(blastMat))
    {
        blastMat[,'length_ratio']=blastMat[,'query_length']/blastMat[,'hit_length']
    }
    #blastMat[,'e.value']=as.numeric(blastMat[,'e.value'])
    subset(blastMat,length_ratio>=threshold)
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

##following function reads protein groups, blast output and filter the protein group and replaces its ids with subject ids in blast result.
proteinGrpFilterReplaceIds<-function(f1,d1,f2,d2,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    blastDB=blastFilterEval(blastFilterMatch(blast(f2,d2),1),upper)
    Mat1=replaceComma(Mat1,3)
    Mat1=upto(Mat1,3,';')
    Mat1=upto(Mat1,3,' ')
    blastDB=upto(blastDB,1,';')
    blastDB=upto(blastDB,1,',')
    blastDB=upto(blastDB,4,' ')
    Mat1[,'protein.accession']=sub("^sw","sp",Mat1[,'protein.accession'])
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    Mat1
}

overlap<-function(replacedIds,Mat)
{
    which(unique(replacedIds) %in% Mat[,'protein.accession'])
}

overlapIndices<-function(replacedIds,Mat)
{
    which(replacedIds %in% Mat[,'protein.accession'])
}


overlap2<-function(replacedIds,Mat)
{
    which(replacedIds %in% Mat)
}

overlapIndices2<-function(replacedIds,Mat)
{
    which(replacedIds %in% Mat)
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

exclusive2<-function(replacedIds,Mat)
{
    which(!replacedIds %in% Mat)
}

voidIsoformEffect<-function(Mat)
{
    Mat=sub("(.*)?\\|([^\\|]+)?\\|(.*)?","\\2",Mat)
    Mat=sub("([^-]+)?-(.*)?","\\1",Mat)
    Mat
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

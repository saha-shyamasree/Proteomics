readPeptides<-function(filename,dir)
{
    read.delim(file=paste(dir,filename,sep=""),header=TRUE)
}
##Filterout peptides that mapps to Reverse and contaminats
peptideFilter<-function(Mat)
{
    Mat[-grep("^CON|^REV",Mat[,'Leading.razor.protein']),]
}

peptideOverlap<-(Mat1,Mat2)
{
    Mat1[which(Mat1[,'Sequence'] %in% Mat2[,'Sequence']),'Sequence']
}

peptideCount<-function(Mat)
{
    length(unique(Mat[,'Sequence'])
}

readProteinGroups<-function(filename,dir)
{
    as.matrix(read.delim(file=paste(dir,filename,sep=""),header=TRUE))
}

proteinGroupCount<-function(Mat)
{
    dim(Mat)  
}

proteinGroupFilter<-function(Mat,threshold)
{
    Mat=Mat[-grep("^CON|^REV",Mat[,'Protein.IDs']),]
    Mat[which(as.numeric(Mat[,'Peptides'])>threshold),]
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

addRazorProteinColumn<-function(Mat)
{
    razorProteins=sub("([^;]+);(.*)?","\\1",Mat[,'Protein.IDs'])
    cbind(Razor.Protein=razorProteins,Mat)
}

##MSGF+ style matrix
newMatrix<-function(Mat)
{
    newMat=matrix(nrow=0,ncol=4)
    colnames(newMat)=c("Protein.Grp","Protein.ID","Peptide.Count","Membership")
    print(dim(Mat)[1])
    for(i in 1:dim(Mat)[1])
    {
        #print(i)
        if(as.numeric(Mat[i,'Number.of.proteins'])>1)
        {
            #print(i)
            proteins=strsplit(Mat[i,'Protein.IDs'],';')
            peptides=strsplit(Mat[i,'Peptide.counts..all.'],';')
            if(length(proteins[[1]])==length(peptides[[1]])) #strsplir returns list
            {
                newMat=rbind(newMat,cbind(rep(i,times=length(proteins[[1]])),proteins[[1]],peptides[[1]],c("anchor",rep("sub",times=length(proteins[[1]])-1))))
            }
            else
            {
                print("ERROR")
            }
        }
        else
        {
            newMat=rbind(newMat,c(i,Mat[i,c('Protein.IDs','Peptide.counts..all.')],"anchor"))
        }
    }
    newMat
}

##replacing any comma in a column value.
replaceComma<-function(Mat,column)
{
    Mat[,column]=gsub(",",";",Mat[,column])
    Mat
}

##may need to manipulate hit_def column to get uniprot id
replaceIds<-function(Mat,B)
{
    ##B should hold blast result in CSV format
    ##Mat is protein identification result from MaxQuant
    
    indices=which(Mat[,'Protein.ID'] %in% B[,'query_name'])
    print(length(indices))
    for(i in indices)
    {
        Mat[i,'Protein.ID'] = as.character(B[which(B[,'query_name'] %in% Mat[i,'Protein.ID']),'hit_def'])
    }
    Mat[,'Protein.ID']
}

filterNewMat<-function(Mat,threshold)
{
    Mat[which(as.numeric(Mat[,'Peptide.Count'])>threshold),]
}
overlap<-function(replacedIds,Mat)
{
    which(unique(replacedIds) %in% Mat[,'Protein.ID'])
}

overlapIndices<-function(replacedIds,Mat)
{
    which(replacedIds %in% Mat[,'Protein.ID'])
}
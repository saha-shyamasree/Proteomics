library(reshape)
library(ggplot2)

readHalfLifeFile<-function(dir,filename)
{
    as.matrix(read.delim(file=paste(dir,filename,sep=""),sep="\t",header=TRUE))
}

drawDensity <-function(Mat)
{
    colnames(mat)=c("uni","trinity")
    #X = rep( 1:dim(Mat)[1], 2)
    Mat_melted = melt(Mat)
    Mat_melted$value = as.numeric( Mat_melted$value)
    #print(colnames(Mat_melted))
    ## Line Plot ##
    p1 <- ggplot(Mat_melted, aes( x = value)) + geom_density(aes(col=factor(X2))) + ggtitle(" Database Length")
    p1
}

prepareDensityMat<-function(MatList,col)
{
    maxLen=0
    for(i in 1:length(MatList))
    {
        maxLen=max(maxLen,dim(MatList[[i]])[1])
    }
    print(maxLen)
    dMat=matrix(nrow=maxLen,ncol=0)
    for(i in 1:length(MatList))
    {
        temp=MatList[[i]][,col]
        print(length(temp))
        length(temp)=maxLen
        print(length(temp))
        dMat=cbind(dMat,temp)
    }
    #drawDensity(dMat)
    dMat
}


#> uniHalfAdeno=readHalfLifeFile("D:/data/Results/Human-Adeno/Half-life/","human_adeno_expasy_uniprot_only.tsv")
#> trinityHalfAdeno=readHalfLifeFile("D:/data/Results/Human-Adeno/Half-life/","trinityPITPRF_compared_uni_half_life_on_seq.tsv")
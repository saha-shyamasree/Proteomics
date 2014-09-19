> f1="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
> f2="only_trinity_PITORF_compare_uni.tsv"
> f3="trinity_PITORF_human_adeno_blast2.csv"
> d1="D:/data/blast/blastCSV/"
> d2="D:/data/Results/Human-Adeno/"
> exc_T=read.table(file=paste(d2,f2,sep=""),header=TRUE,sep="\t")
> Uni=read.csv(file=paste(d2,f1,sep=""),header=TRUE)
> blastRes=blast(f3,d1)



combineRes<-function(exc,Mat,blastRes)
{
    exc=replaceComma(exc,3)
    exc=replaceComma(exc,5)
    Mat[,'protein.accession']=sub("^sw","sp",Mat[,'protein.accession'])
    Mat=replaceComma(Mat,3)
    Mat=replaceComma(Mat,5)
    combinedMat=matrix(nrow=0,ncol=(dim(exc)[2]+dim(Mat)[2]+dim(blastRes)[2]),dimnames=list(c(),c(colnames(exc),colnames(Mat),colnames(blastRes))))
    for(i in 1:dim(exc)[1])
    {
        print(paste("i:",i))
        ind1=which(Mat[,'protein.accession']==exc[i,'protein.accession'] )
        if(length(ind1)>0)
        {
            print(paste("ind1:",ind1))
            ind2=which(blastRes[,'query_name']==exc[i,'description'])
            if(length(ind2)>0)
            {
                print(paste("ind2:",ind2))
                tcol1=cbind(exc[i,],Mat[ind1,],deparse.level=0)
                tcol=cbind(tcol1,blastRes[ind2,],deparse.level=0)
                combinedMat=rbind(combinedMat,tcol)
            }
            else
            {
                print("ERROR")
            }
        }
        else
        {
            ind2=which(blastRes[,'query_name']==exc[i,'description'])
            if(length(ind2)>0)
            {
                print(paste("ind2:",ind2))
                empty=t(matrix(rep.int(0,dim(Mat)[2])))
                colnames(empty)=colnames(Mat)
                tcol1=cbind(exc[i,],empty,deparse.level=0)
                tcol=cbind(tcol1,blastRes[ind2,],deparse.level=0)
                combinedMat=rbind(combinedMat,tcol)
            }
            else
            {
                print("ERROR")
            }
        }
    }
    write.csv(combinedMat,file=paste(d2,"combinedTrinityOnlyWithBlastAndMatprot.csv"),quote = FALSE,row.names = FALSE)
}

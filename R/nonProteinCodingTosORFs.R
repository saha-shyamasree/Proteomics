## This R code check how many of non-protein coding ORFs are sORFs

#### Common values
upper=0.000000000000000000000000000001
pepThreshold=1
rev=1
peptide=1
len=66
################## PASA non-protein coding class ##################
f1="nonCodingPrt.tsv"
d1="D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/"
f3="human_adeno_mydb_pasa.assemblies_ORFsV1.csv"
d3="D:/data/blast/blastCSV/PASA/Human-Adeno/"
Mats=readMatrices(f1,f3,d1,d3,rev,peptide,pepThreshold)
sORFCountSingleDB(Mats$Mat1, Mats$blastDB1, Mat1Name, upper, len)


proteinGroup <- function(filename,dirc,sep)
{
    as.matrix(read.table(file=paste(dirc,filename,sep=""),sep=sep,header=TRUE))
}
readMatrices<-function(f1,f3,d1,d3,rev=1,peptide=1,pepThreshold=1)
{
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1,'\t'),rev,peptide,pepThreshold)
    
    blastDB1=blast(f3,d3)

    Mat1=replaceComma(Mat1,3)
    Mat1=upto(Mat1,3,';')
    Mat1=upto(Mat1,3,' ')
    blastDB1=upto(blastDB1,1,';')
    blastDB1=upto(blastDB1,1,',')
    blastDB1=upto(blastDB1,4,' ')
    
    Mat1[,'protein.accession']=sub("^sw","sp",Mat1[,'protein.accession'])

    list(Mat1=Mat1, blastDB1=blastDB1)
}

sORFCountSingleDB<-function(Mat1, blastDB1, Mat1Name, upper, len)
{
    Mat1IdentifiedBlast=blastDB1[which(blastDB1[,'query_name'] %in% Mat1[,'protein.accession']),]
    
    sORFMat1=Mat1IdentifiedBlast[which(Mat1IdentifiedBlast[,'query_length']<len),]
    
    print("Mat1 sORF")
    print(dim(sORFMat1))
    
    print("Mat1 sORF: 3'uORF")
    print(length(grep("^Dataset_C",sORFMat1[,'query_name'])))
    
    print("Mat1 sORF: 5'uORF")
    print(length(grep("^Dataset_B",sORFMat1[,'query_name'])))

}
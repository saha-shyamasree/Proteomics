#### Conrad's proteomics course.

### This code reads identified ORFs and BLAST result and write a list of uniprot proteins that were identified.
source("D:/Code/Proteomics/R/RLib.R")
rev=1
peptide=1
pepThreshold=1
upper=0.000000000000000000000000000001

f1="trinity_PITORF+fdr+th+grouping+prt.csv"
f2="trinity_PITORF_human_adeno_blast2.csv"
d1="D:/data/Results/Human-Adeno/GIOPaperResults/trinity/"
d2="D:/data/blast/blastCSV/"

readMat<- function(f1,f2,d1,d2,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    #print(f1)
    #print(f2)
    #print(d1)
    #print(d2)
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    Mat1=replaceComma(Mat1,3)
    Mat1[,'protein.accession']=sub("^sw","sp",Mat1[,'protein.accession'])
    
    blastDB=blastFilterEval(blastFilterMatch(blast(f2,d2),1),upper)
    print("Mat1")
    print(dim(Mat1))
    print("BLAST")
    print(dim(blastDB))
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    Mat1[,'protein.accession']=sub("^sp","sw",Mat1[,'protein.accession'])
    list('Mat1'=Mat1,'blast'=blastDB)
}


write.table(Mats$Mat1[,'protein.accession'],file=paste(d1,"IdentifiedORFsUniprotIds.tsv"),col.names=F,row.names=F,quote=F)

## This R code try to find relation between non-protein coding peptides identified in a PIT approaced search against a standard
## database search.

source("D:/Code/Proteomics/R/RLib.R")
readFile<-function(filename,dir,sep)
{
    as.matrix(read.table(file=paste(dir,filename,sep=""),sep=sep, header=TRUE))
}

fname1="DM_from_raw_uniprot+fdr+th+grouping.csv"
fname1Prt="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
dir1="D:/data/Results/Human-Adeno/Identification/"

fname2="nonCodingPSM.tsv"
fname2Prt="nonCodingPrt.tsv"
dir2="D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/"

stdPrt=readFile(fname1Prt,dir1,',')
stdPSM=readFile(fname1,dir1,',')
stdPep=unique(stdPSM[,'Sequence'])

nonCDSPrt=readFile(fname2Prt,dir2,'\t')
nonCDSPSM=readFile(fname2,dir2,'\t')
nonCDSPep=unique(nonCDSPSM[,'Sequence'])

overlapPep=nonCDSPep[which(nonCDSPep %in% stdpep)]

stdProtOverlap=unique(stdPSM[which(stdPSM[,'Sequence'] %in% overlap),'proteinacc_start_stop_pre_post_.'])
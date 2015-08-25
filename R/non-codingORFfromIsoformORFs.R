## This R code Reads (1) the list of known proteins with/without SAPs, (2) isoforms with/without SAPs, and (3) the list of ORFs classified according to their
## transcript bio-type. Aim of this code is to exclude (1) from the list of (3).

source("D:/Code/Proteomics/R/RLib.R")
readList<-function(filepath)
{
    as.matrix(read.csv(file=filepath, header=TRUE))
}


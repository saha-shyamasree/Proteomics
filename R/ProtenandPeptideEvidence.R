
peptideFile="trinity_PITORF+fdr+th+grouping.csv"
peptideFolder="D:/data/Results/Human-Adeno/Identification/"
proteinFile="TrinityOnly_novelProteinsCandidate.csv"
proteinFolder="D:/data/Results/Human-Adeno/Identification/"

readPeptideFile <- function(filename, directory)
{
    as.matrix(read.csv(file=paste(directory,filename,sep=""),header=TRUE))
}

readProteinFile <- function(filename, directory)
{
    as.matrix(read.csv(file=paste(directory,filename,sep=""),header=TRUE))
}

peptideSearch<-function(proteinId, peptideMat, col)
{
    grep(proteinId,peptideMat[,col])
}

mergeProteinAndPeptideEvidences<-function(proteinFile, proteinDir, peptideFile, peptideDir)
{
    proteins=readProteinFile(proteinFile, proteinDir)
    peptides=readPeptideFile(peptideFile, peptideDir)
    mergedEvidences=NULL
    for(i in 1:dim(proteins)[1])
    {
        indxs=peptideSearch(proteins[i,'protein.accession'],peptides,20)
        mergedEvidences=rbind(mergedEvidences,cbind(do.call('rbind',replicate(length(indxs),proteins[i,],simplify=FALSE)),peptides[indxs,]))
    }
    mergedEvidences
}

 mergeProteinAndPeptideEvidences(proteinFile,proteinFolder,peptideFile,peptideFolder)
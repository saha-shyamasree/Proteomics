##This code takes list of known proteins, known proteins with SAPs and possible isoforms produced by 'IdentifyProteinIsoformSAP.py' and
##identify how many of these are found from peptide and protein identification.

source("D:/Code/Proteomics/R/RLib.R")
readList<-function(filepath)
{
    as.matrix(read.csv(file=filepath, header=TRUE))
}

#identifiedORFs is filtered and protein accession replaced by uniprot ids matrix of identified ORFs, knownProteinList is list of ORFs with exact match to an Uniprot protein
identifiedORFsClassKnownProtein<-function(identifiedORFs, knownProteinList, outFile)
{
    Mat=identifiedORFs[which(identifiedORFs[,'description'] %in% knownProteinList[,'ORF.Id']),]
    write.table(Mat,file=outFile,sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
}

identifiedORFsClassKnownProteinSAP<-function(identifiedORFs, knownProteinSAPList, outFile)
{
    Mat=identifiedORFs[which(identifiedORFs[,'description'] %in% knownProteinSAPList[,'ORF.Id']),]
    write.table(Mat,file=outFile,sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
}

identifiedORFsClassIsoform<-function(identifiedORFs, isoformList, outFile)
{
    Mat=identifiedORFs[which(identifiedORFs[,'description'] %in% isoformList[,'ORF.Id']),]
    index=match(Mat[,'description'],isoformList[,'ORF.Id'])
    Mat=cbind(Mat,isoformList[index,'Type'])
    write.table(Mat,file=outFile,sep='\t',quote = FALSE,row.names = FALSE, col.names=TRUE)
}



################################# 20-05-2015 ##############################################
#### PASA assembly for human adenovirus

##File Paths and other thresholding values
identifiedORFsDir="D:/data/Results/Human-Adeno/Identification/PASA/sORF/"
identifiedORFsFilename="pasa_assemblyV1+fdr+th+grouping+prt.csv"

knownProteinPath="D:/data/blast/blastCSV/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsV2.csv"
knownProteinSAPPath="D:/data/blast/blastCSV/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsSAPsV2.csv"
isoformPath="D:/data/blast/blastCSV/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsV2.csv"
isoformSAPsPath="D:/data/blast/blastCSV/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsSAPsV2.csv"

knownIdentifiedProteinsPath=paste(identifiedORFsDir,"pasa_assemblyV1_knownProteinsV2.tsv")
knownSAPIdentifiedProteinsPath=paste(identifiedORFsDir,"pasa_assemblyV1_knownSAPProteinsV2.tsv")
isoformIdentifiedProteinsPath=paste(identifiedORFsDir,"pasa_assemblyV1_isoformProteinsV2.tsv")
isoformSAPsIdentifiedProteinsPath=paste(identifiedORFsDir,"pasa_assemblyV1_isoformSAPsProteinsV2.tsv")

blastFile="human_adeno_mydb_pasa.assemblies_ORFsV1.csv"
blastDir="D:/data/blast/blastCSV/"
upper=0.000000000000000000000000000001

##Reading matrices
identifiedORFs=proteinGrpFilterReplaceIds(identifiedORFsFilename,identifiedORFsDir,blastFile,blastDir,1,1,1,upper)
knownProteinList=readList(knownProteinPath)
knownProteinSAPList=readList(knownProteinSAPPath)
isoformList=readList(isoformPath)
isoformSAPsList=readList(isoformSAPsPath)

identifiedORFsClassKnownProtein(identifiedORFs,knownProteinList,knownIdentifiedProteinsPath)
identifiedORFsClassKnownProteinSAP(identifiedORFs,knownProteinSAPList,knownSAPIdentifiedProteinsPath)
identifiedORFsClassIsoform(identifiedORFs,isoformList,isoformIdentifiedProteinsPath)
identifiedORFsClassIsoform(identifiedORFs,isoformSAPsList,isoformSAPsIdentifiedProteinsPath)
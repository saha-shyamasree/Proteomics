readPSM<-function(filename,dir)
{
    as.matrix(read.csv(file=paste(dir,filename,sep=""),header=TRUE))
}

countPSMs<-function(Mat1)
{
    length(unique(Mat1[,'Spectrum.ID']))
}

overlapPepIndices<-function(Mat1,Mat2)
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
filename1="trinity_PITORF+fdr+th+grouping.csv"
filename2="DM_from_raw_uniprot+fdr+th+grouping.csv"
folder="C:/Users/shyama/Dropbox/results/Human_adenovirus/"
TAdenoPSM=readPSM(filename1,folder)
print("trinity mat dim and peptide count")
dim(TAdenoPSM)
length(unique(TMat[,11]))
UAdenoPSM=readPSM(filename2,folder)
print("uniprot mat dim and peptide count")
dim(UAdenoPSM)
length(unique(UMat[,11]))
print("Overlap")
length(overlap(TMat,UMat,11,11))

#################### Human Hendra ########################
filenameHHT="trinity+fdr+th+grouping.csv"
filenameHHU="uniprot_human+fdr+th+grouping.csv"
folder="C:/Users/shyama/Dropbox/results/Bat_human_hendra/Human/"
THHPSM=readPSM(filenameHHT,folder)
UHHPSM=readPSM(filenameHHU,folder)
length(unique(THHPSM[,'Sequence']))
length(unique(UHHPSM[,'Sequence']))
length(unique(THHPSM[,'Spectrum.ID']))
length(unique(UHHPSM[,'Spectrum.ID']))
length(which(unique(THHPSM[,'Spectrum.ID']) %in% unique(UHHPSM[,'Spectrum.ID'])))
length(which(unique(THHPSM[,'Sequence']) %in% unique(UHHPSM[,'Sequence'])))

##################### Bat Hendra ##########################
filenameBHT="trinity_bat+fdr+th+grouping.csv"
filenameBHU="uniprot_bat+fdr+th+grouping.csv"
folder="C:/Users/shyama/Dropbox/results/Bat_human_hendra/Bat/"
TBHPSM=readPSM(filenameBHT,folder)
UBHPSM=readPSM(filenameBHU,folder)
length(unique(TBHPSM[,'Sequence']))
length(unique(UBHPSM[,'Sequence']))
length(unique(TBHPSM[,'Spectrum.ID']))
length(unique(UBHPSM[,'Spectrum.ID']))
length(which(unique(TBHPSM[,'Spectrum.ID']) %in% unique(UBHPSM[,'Spectrum.ID'])))
length(which(unique(TBHPSM[,'Sequence']) %in% unique(UBHPSM[,'Sequence'])))

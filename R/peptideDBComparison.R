##This code reads two tsv files, two peptide identification tsv from MSGF+ and mzidenML-lib. Use case, peptides from Standard and PIT searches
##and peptide DBs for these searches [generated using Jun Fan's digest code and and those files were converted from fasta to tabular format].

source("D:/Code/Proteomics/R/RLib.R")

readPeptides<-function(filename,dir)
{
    as.matrix(read.table(file=paste(dir,filename,sep=""), sep="\t", header=FALSE))
}

overlapPepVector<-function(Mat1,Mat2)
{
    which(Mat1 %in% Mat2)
}

nonoverlapPepVector<-function(Mat1,Mat2)
{
    which(! Mat1 %in% Mat2)
}

readFiles<-function(f1, f2, f3, f4, d1, d2, d3, d4)
{
    Mat1PSM=readPSM(f1,d1)
    Mat2PSM=readPSM(f2,d2)
    
    Mat1Peptides=readPeptides(f3,d3)
    Mat2Peptides=readPeptides(f4,d4)
    
    list("Mat1PSM"=Mat1PSM,"Mat2PSM"=Mat1PSM,"Mat1Peptides"=Mat1Peptides,"Mat2Peptides"=Mat2Peptides)
}

checkIdentifiedPeptidesInPeptideDB <- function(Mat1PSM, Mat1Peptides)
{
    Mat1PSMUnq=unique(Mat1PSM[,column1])
    Mat1PSMUnq[which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2]))]
    Mat1PSMUnqMetAdded=paste("M",Mat1PSMUnq[which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2]))],sep="")
    Mat1PSMUnq[which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2]))]=Mat1PSMUnqMetAdded
    length(which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2])))##this should return 0.
}

checkExistenceofSearchSpecificPeptidesinOtherDB <- function(Mat1PSM, Mat2Peptides)
{
    
}

dir1="D:/data/Results/Human-Adeno/GIOPaperResults/trinity/"
dir2="D:/data/Results/Human-Adeno/GIOPaperResults/Uniprot-trinity/"
dir3="D:/data/Data/fasta/GIOPaper/"
dir4=dir3
fname1="trinity_PITORF+fdr+th+grouping.csv"
fname2="DM_from_raw_uniprot+fdr+th+grouping.csv"
fname3="Supplementary_dataset_1_-_human_trinity_assembled_transcripts-ORF_concatenated_target_decoyPeptidesv4.tsv"
fname4="human_adenovirus_concatenated_target_decoyPeptidesv4.tsv"
column1=11
column2=2
#peptideDBCheck(fname1, dir1, fname2, dir2, fname3, dir3, fname4, dir4, column1, column2)

#peptideDBCheck<-function(fname1, dir1, fname2, dir2, fname3, dir3, fname4, dir4, column1, column2)
#{
Mat1PSM=readPSM(fname1,dir1)
Mat2PSM=readPSM(fname2,dir2)

Mat1Peptides=readPeptides(fname3,dir3)
Mat2Peptides=readPeptides(fname4,dir4)

Mat1PSMUnq=unique(Mat1PSM[,column1])
Mat1PSMUnq[which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2]))]
Mat1PSMUnqMetAdded=paste("M",Mat1PSMUnq[which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2]))],sep="")
Mat1PSMUnq[which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2]))]=Mat1PSMUnqMetAdded
length(which(! Mat1PSMUnq %in% unique(Mat1Peptides[,column2])))##this should return 0.

length(which(! Mat1PSMUnq %in% unique(Mat2Peptides[,column2])))
length(which(! Mat2PSMUnq %in% unique(Mat1Peptides[,column2])))

Mat2PSMUnq=unique(Mat2PSM[,column1])
Mat2PSMUnq[which(! Mat2PSMUnq %in% unique(Mat2Peptides[,column2]))]

##This is to check that Jun's code has generated at least the peptides that have been identified by the MSGF+ tool.
M1Len=length(overlapPepVector(Mat1PSMUnq,Mat1Peptides[,column2]))
M2Len=length(overlapPepVector(Mat2PSMUnq,Mat2Peptides[,column2]))

print(head(unique(Mat1PSM[,column1])))
print(head(unique(Mat2PSM[,column1])))
print(head(Mat1Peptides))
print(head(Mat2Peptides))

if(M1Len==length(unique(Mat1PSM[,column1])))
{
    print("All identified Peptides Found in Mat1 are in the Peptide list")
    if(M2Len==length(unique(Mat2PSM[,column1])))
    {
        print("All identified Peptides Found in Mat2 are in the Peptide list")
    }
    else
    {
        print("All identified Peptides are not found in Mat2 Peptide list")
        print(M2Len)
        print(length(unique(Mat2PSM[,column1])))
        print(length(unique(Mat2Peptides[,column2])))
    }
}else
{
    print("All identified Peptides are not found in Mat1 peptide list")
    print(M1Len)
    print(length(unique(Mat1PSM[,column1])))
    print(length(unique(Mat1Peptides[,column2])))
}
#}

################################ Human Adeno virus standard vs genome guided PIT (cufflinks) ######################################

dirc="D:/data/Results/Human-Adeno/GIOPaperResults/cufflinks/"
diru="D:/data/Results/Human-Adeno/GIOPaperResults/Uniprot-Cufflinks/"
dircomm="D:/data/Data/fasta/GIOPaper/"

fnamec="cufflinks_main_PITORF+fdr+th+grouping.csv"
fnameu="human_uniprot+fdr+th+grouping.csv"
fnamec2="cufflinks_main-ORF_concatenated_target_decoyPeptidesv4.tsv"
fnameu2="HUMAN_concatenated_target_decoyPeptidesv4.tsv"
column1=11
column2=2

MatCPSM=readPSM(fnamec,dirc)
MatUPSM=readPSM(fnameu,diru)

MatCPeptides=readPeptides(fnamec2,dircomm)
MatUPeptides=readPeptides(fnameu2,dircomm)

MatCPSMUnq=unique(MatCPSM[,column1])
head(MatCPSMUnq[which(! MatCPSMUnq %in% unique(MatCPeptides[,column2]))])
MatCPSMUnqMetAdded=paste("M",MatCPSMUnq[which(! MatCPSMUnq %in% unique(MatCPeptides[,column2]))],sep="")
MatCPSMUnq[which(! MatCPSMUnq %in% unique(MatCPeptides[,column2]))]=MatCPSMUnqMetAdded
length(which(! MatCPSMUnq %in% unique(MatCPeptides[,column2])))##this should return 0.

MatUPSMUnq=unique(MatUPSM[,column1])
head(MatUPSMUnq[which(! MatUPSMUnq %in% unique(MatUPeptides[,column2]))])
MatUPSMUnqMetAdded=paste("M",MatUPSMUnq[which(! MatUPSMUnq %in% unique(MatUPeptides[,column2]))],sep="")
MatUPSMUnq[which(! MatUPSMUnq %in% unique(MatUPeptides[,column2]))]=MatUPSMUnqMetAdded
length(which(! MatUPSMUnq %in% unique(MatUPeptides[,column2])))##this should return 0.


##This is to check that Jun's code has generated at least the peptides that have been identified by the MSGF+ tool.
MCLen=length(overlapPepVector(MatCPSMUnq,MatCPeptides[,column2]))
MULen=length(overlapPepVector(MatUPSMUnq,MatUPeptides[,column2]))

length(which(! MatCPSMUnq %in% unique(MatUPeptides[,column2])))
length(which(! MatUPSMUnq %in% unique(MatCPeptides[,column2])))

MatCPSMUnq[which(! MatCPSMUnq %in% unique(MatUPeptides[,column2]))]
MatUPSMUnq[which(! MatUPSMUnq %in% unique(MatCPeptides[,column2]))]

### count of Uniquely identified peptides that are not in the opposite database
length(which(! Mat1PepUniqSeq %in% unique(Mat2Peptides[,column2]))) ##MatUPepUniqSeq from OverlappingPeptideVsProteins.R
length(which(! Mat2PepUniqSeq %in% unique(Mat1Peptides[,column2]))) ##MatCPepUniqSeq from OverlappingPeptideVsProteins.R
length(which(! MatUPepUniqSeq %in% unique(MatCPeptides[,column2]))) ##MatUPepUniqSeq from OverlappingPeptideVsProteins.R
length(which(! MatCPepUniqSeq %in% unique(MatUPeptides[,column2]))) ##MatCPepUniqSeq from OverlappingPeptideVsProteins.R

### Peptide list of Uniquely identified peptides that are in the opposite database
Mat1PepUniqSeqDBIn = unique(Mat1PepUniqSeq[which( Mat1PepUniqSeq %in% unique(Mat2Peptides[,column2]))])
Mat2PepUniqSeqDBIn = unique(Mat2PepUniqSeq[which( Mat2PepUniqSeq %in% unique(Mat1Peptides[,column2]))])
MatUPepUniqSeqDBIn = unique(MatUPepUniqSeq[which( MatUPepUniqSeq %in% unique(MatCPeptides[,column2]))])
MatCPepUniqSeqDBIn = unique(MatCPepUniqSeq[which( MatCPepUniqSeq %in% unique(MatUPeptides[,column2]))])


### Uniquely identified peptides(even though they exist in the opposite DB) contributing to one-hit
length(which( Mat1PepUniqSeqDBIn %in% Mat1PrtSinglePepUniq)) ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R
length(which( Mat2PepUniqSeqDBIn %in% Mat2PrtSinglePepUniq)) ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R
length(which( MatCPepUniqSeqDBIn %in% MatCPrtSinglePepUniq)) ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R
length(which( MatUPepUniqSeqDBIn %in% MatUPrtSinglePepUniq)) ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R

### Uniquely identified peptides (even though they exists in the other DB) that are not in one hit one wonder.
Mat1PepUniqSeqDBInBothDB=Mat1PepUniqSeqDBIn[ which(! Mat1PepUniqSeqDBIn %in% Mat1PrtSinglePepUniq)] ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R
Mat2PepUniqSeqDBInBothDB=Mat2PepUniqSeqDBIn[ which(! Mat2PepUniqSeqDBIn %in% Mat2PrtSinglePepUniq)] ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R
MatCPepUniqSeqDBInBothDB=MatCPepUniqSeqDBIn[ which(! MatCPepUniqSeqDBIn %in% MatCPrtSinglePepUniq)] ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R
MatUPepUniqSeqDBInBothDB=MatUPepUniqSeqDBIn[ which(! MatUPepUniqSeqDBIn %in% MatUPrtSinglePepUniq)] ##MatCPrtSinglePepUniq from OverlappingPeptideVsProteins.R

############################## Excluding One hit wonder ##########################################
### Check the PSMs mapping to the peptides that exist in both database but identified only in one.

Mat1PSMInBothDB=Mat1PSM[which(Mat1PSM[,'Sequence'] %in% Mat1PepUniqSeqDBInBothDB),]
Mat2PSMInBothDB=Mat2PSM[which(Mat2PSM[,'Sequence'] %in% Mat2PepUniqSeqDBInBothDB),]
MatCPSMInBothDB=MatCPSM[which(MatCPSM[,'Sequence'] %in% MatCPepUniqSeqDBInBothDB),]
MatUPSMInBothDB=MatUPSM[which(MatUPSM[,'Sequence'] %in% MatUPepUniqSeqDBInBothDB),]

### Check the peptides in the other DB supported by the above PSMs.

Mat1PSMDiffMapInOtherDB=Mat1PSM[which(Mat1PSM[,'Spectrum.ID'] %in% Mat2PSMInBothDB[,'Spectrum.ID']),]
Mat2PSMDiffMapInOtherDB=Mat2PSM[which(Mat2PSM[,'Spectrum.ID'] %in% Mat1PSMInBothDB[,'Spectrum.ID']),]
MatCPSMDiffMapInOtherDB=MatCPSM[which(MatCPSM[,'Spectrum.ID'] %in% MatUPSMInBothDB[,'Spectrum.ID']),]
MatUPSMDiffMapInOtherDB=MatUPSM[which(MatUPSM[,'Spectrum.ID'] %in% MatCPSMInBothDB[,'Spectrum.ID']),]



length(unique(Mat1PSMDiffMapInOtherDB[,'Spectrum.ID']))
length(unique(Mat2PSMDiffMapInOtherDB[,'Spectrum.ID']))
length(unique(MatCPSMDiffMapInOtherDB[,'Spectrum.ID']))
length(unique(MatUPSMDiffMapInOtherDB[,'Spectrum.ID']))

length(unique(Mat1PSMDiffMapInOtherDB[,'Sequence']))
length(unique(Mat2PSMDiffMapInOtherDB[,'Sequence']))
length(unique(MatCPSMDiffMapInOtherDB[,'Sequence']))
length(unique(MatUPSMDiffMapInOtherDB[,'Sequence']))

######################################## All uniquely identified peptides that exist in the other DB but was not identified for that DBB search
### Check the PSMs mapping to the peptides that exist in both database but identified only in one.

Mat1PSMInBothDB=Mat1PSM[which(Mat1PSM[,'Sequence'] %in% Mat1PepUniqSeqDBIn),]
Mat2PSMInBothDB=Mat2PSM[which(Mat2PSM[,'Sequence'] %in% Mat2PepUniqSeqDBIn),]
MatCPSMInBothDB=MatCPSM[which(MatCPSM[,'Sequence'] %in% MatCPepUniqSeqDBIn),]
MatUPSMInBothDB=MatUPSM[which(MatUPSM[,'Sequence'] %in% MatUPepUniqSeqDBIn),]

length(unique(Mat1PSMInBothDB[,'Spectrum.ID']))
length(unique(Mat2PSMInBothDB[,'Spectrum.ID']))
length(unique(MatCPSMInBothDB[,'Spectrum.ID']))
length(unique(MatUPSMInBothDB[,'Spectrum.ID']))

### Check the peptides in the other DB supported by the above PSMs.

Mat1PSMDiffMapInOtherDB=Mat1PSM[which(Mat1PSM[,'Spectrum.ID'] %in% Mat2PSMInBothDB[,'Spectrum.ID']),]
Mat2PSMDiffMapInOtherDB=Mat2PSM[which(Mat2PSM[,'Spectrum.ID'] %in% Mat1PSMInBothDB[,'Spectrum.ID']),]
MatCPSMDiffMapInOtherDB=MatCPSM[which(MatCPSM[,'Spectrum.ID'] %in% MatUPSMInBothDB[,'Spectrum.ID']),]
MatUPSMDiffMapInOtherDB=MatUPSM[which(MatUPSM[,'Spectrum.ID'] %in% MatCPSMInBothDB[,'Spectrum.ID']),]



length(unique(Mat1PSMDiffMapInOtherDB[,'Spectrum.ID']))
length(unique(Mat2PSMDiffMapInOtherDB[,'Spectrum.ID']))
length(unique(MatCPSMDiffMapInOtherDB[,'Spectrum.ID']))
length(unique(MatUPSMDiffMapInOtherDB[,'Spectrum.ID']))

length(unique(Mat1PSMDiffMapInOtherDB[,'Sequence']))
length(unique(Mat2PSMDiffMapInOtherDB[,'Sequence']))
length(unique(MatCPSMDiffMapInOtherDB[,'Sequence']))
length(unique(MatUPSMDiffMapInOtherDB[,'Sequence']))

##This code compares identified peptides from same MS data against different protein databases and check whether same peptide is in different sides of the threshold in different databse search.

source("D:/Code/Proteomics/R/RLib.R")
### GIO Paper Analysis ###

### De-novo vs standard
tFile="trinity_PITORF+fdr+th+grouping2.csv"
tDir="D:/data/Results/Human-Adeno/GIOPaperResults/trinity/"
uFile="DM_from_raw_uniprot+fdr+th+grouping2.csv"
uDir="D:/data/Results/Human-Adeno/GIOPaperResults/Uniprot-trinity/"

##Reading PSM files
tPSMs=readPSM(tFile, tDir)
uPSMs=readPSM(uFile, uDir)

##Separating PSMs that passed the threshold
tPSMPassed=thresholdPassedPSMs(tPSMs)
uPSMPassed=thresholdPassedPSMs(uPSMs)

##Separating PSMs that did not passe the threshold
tPSMFailed=belowThresholdPSMs(tPSMs)
uPSMFailed=belowThresholdPSMs(uPSMs)

##Unique peptide sequences identified by passed PSMs
tPSMPassedUniqueSeq=unique(tPSMPassed[,'Sequence'])
uPSMPassedUniqueSeq=unique(uPSMPassed[,'Sequence'])

##Unique peptide sequences identified by failed PSMs
tPSMFailedUniqueSeq=unique(tPSMFailed[,'Sequence'])
uPSMFailedUniqueSeq=unique(uPSMFailed[,'Sequence'])

##Check whether any passed trinity sequence overlaps with failed uniprot sequence
tFailedinUPassed=tPSMFailedUniqueSeq[which(tPSMFailedUniqueSeq %in% uPSMPassedUniqueSeq)]

##Check whether any passed uniprot sequence overlaps with failed trinity sequence
uFailedinTPassed=uPSMFailedUniqueSeq[which(uPSMFailedUniqueSeq %in% tPSMPassedUniqueSeq)]

### Genome-Guided vs standard
cFile="cufflinks_main_PITORF+fdr+th+grouping2.csv"
cDir="D:/data/Results/Human-Adeno/GIOPaperResults/cufflinks/"
uFile="human_uniprot+fdr+th+grouping2.csv"
uDir="D:/data/Results/Human-Adeno/GIOPaperResults/Uniprot-Cufflinks/"

##Reading PSM files
cPSMs=readPSM(cFile, cDir)
uPSMs=readPSM(uFile, uDir)

##Separating PSMs that passed the threshold
cPSMPassed=thresholdPassedPSMs(cPSMs)
uPSMPassed=thresholdPassedPSMs(uPSMs)

##Separating PSMs that did not passe the threshold
cPSMFailed=belowThresholdPSMs(cPSMs)
uPSMFailed=belowThresholdPSMs(uPSMs)

##Unique peptide sequences identified by passed PSMs
cPSMPassedUniqueSeq=unique(cPSMPassed[,'Sequence'])
uPSMPassedUniqueSeq=unique(uPSMPassed[,'Sequence'])

##Unique peptide sequences identified by failed PSMs
cPSMFailedUniqueSeq=unique(cPSMFailed[,'Sequence'])
uPSMFailedUniqueSeq=unique(uPSMFailed[,'Sequence'])

##Check whether any passed cufflinks sequence overlaps with failed uniprot sequence
cFailedinUPassed=cPSMFailedUniqueSeq[which(cPSMFailedUniqueSeq %in% uPSMPassedUniqueSeq)]
length(cFailedinUPassed)

##Check whether any passed uniprot sequence overlaps with failed cufflinks sequence
uFailedinCPassed=uPSMFailedUniqueSeq[which(uPSMFailedUniqueSeq %in% cPSMPassedUniqueSeq)]
length(uFailedinCPassed)
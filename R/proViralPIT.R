### This code uses protein overlap code and count the host protein and the virus protein identified from the PIT search.

source("/data/home/btw796/Code2/Proteomics/R/RLib.R")
#sink(file="/data/home/btw796/Code2/Proteomics/R/proViralLog.txt")

replaceIds<-function(Mat,B)
{
    ##B should hold blast result in CSV format
    ##Mat is protein identification result from MSGF+
    indices=which(Mat[,'description'] %in% B[,'query_name'])
    print("blast map count")
    print(length(indices))
    for(i in indices)
    {
        Mat[i,'description'] = as.character(B[which(B[,'query_name'] %in% Mat[i,'description']),'hit_def'])
    }
    Mat[,'description']
}


identifiedBlastResult<-function(Mat,B,out)
{
    ##B should hold blast result in CSV format
    ##Mat is protein identification result from MSGF+
    BIdentified=B[which(B[,'query_name'] %in% Mat[,'description']),]
    write.csv(BIdentified,file=out, quote = FALSE,row.names = FALSE)
}

readMats<- function(f1,f2,d1,d2,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    #print(f1)
    #print(f2)
    #print(d1)
    #print(d2)
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    blastDB=blastFilterEval(blastFilterMatch(blast(f2,d2),1),upper)
    print("Mat1")
    print(dim(Mat1))
    print("BLAST")
    print(dim(blastDB))
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    
    write.csv(Mat1,file=paste(d1,f1,"DBIds.csv",sep=""))
    identifiedBlastResult(Mat1, blastDB, paste(d2,"Identified/",f2,sep=""))
}

readMat<- function(f1,f2,d1,d2,rev=1,peptide=1,pepThreshold=1,upper=0.1)
{
    #print(f1)
    #print(f2)
    #print(d1)
    #print(d2)
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold)
    blastDB=blastFilterEval(blastFilterMatch(blast(f2,d2),1),upper)
    print("Mat1")
    print(dim(Mat1))
    print("BLAST")
    print(dim(blastDB))
    print("No of ORFs:")
    print(dim(blast(f2,d2)))
    Mat1=replaceComma(Mat1,5)
    #print(head(Mat1[,'protein.accession']))
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    #print(head(Mat1[,'protein.accession']))
    list('Mat1'=Mat1,'blast'=blastDB)
}

PAGCount<-function(Mat)
{
	seenProt=c()
	seenPAG=c()
	count=0
	bMat=Mat
	for(i in 1:dim(Mat)[1])
	{
		if(Mat[i,'protein.accession'] %in% seenProt)
		{
			if(Mat[i,'PAG.ID'] %in% seenPAG)
			{}
			else
			{
				seenPAG=c(seenPAG,Mat[i,'PAG.ID'])
				
			}
						
		}
		else
		{
			if(Mat[i,'PAG.ID'] %in% seenPAG)
			{
				seenProt=c(seenProt,Mat[i,'protein.accession'])
			}
			else
			{
				seenProt=c(seenProt,Mat[i,'protein.accession'])
				seenPAG=c(seenPAG,Mat[i,'PAG.ID'])
				count=count+1
			}
		}	
	}
	print("PAG Count")
	print(count)
}
rev=1
peptide=1
pepThreshold=1
upper=0.000000000000000000000000000001
revStr="_REVERSED"
contStr=""

#########  Bat Hendra ###################
f1="trinity_bat+fdr+th+grouping+prt.csv"
f2="Trinity_P_alecto-ORF.fasta_reviewed.xml2.csv"
#f2="Trinity_P_alecto-ORF.fasta_All.xml2.csv"
f3="trinity_bat+fdr+th+grouping.csv"
dir1="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/Hendra/identification/"
dir2="/data/SBCS-BessantLab/shyama/Data/Bristol/Signature_peptide_seeded_database/Bat/"
host="Pteropus alecto"
virus="_HENDH"

########## Human Adenovirus ###################

#f1="trinity_PITORF+fdr+th+grouping+prt.csv"
#f2="Supplementary_dataset_1_-_human_trinity_assembled_transcripts-ORF.fasta_reviewed.xml2.csv"
#f2="Supplementary_dataset_1_-_human_trinity_assembled_transcripts-ORF.fasta_All.xml2.csv"
#f3="trinity_PITORF+fdr+th+grouping.csv"
#dir1="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/"
#dir2="/data/SBCS-BessantLab/shyama/Data/Bristol/Signature_peptide_seeded_database/Human/"
#host="_HUMAN"
#virus="_ADE05"


############## Mosquito ######################

#f1="mosquito+fdr+th+grouping+prt.csv"
#f2="aedesTrinity_ORF.fasta_reviewed.xml2.csv"
#f2="aedesTrinity_ORF.fasta_All.xml2.csv"
#f3="mosquito+fdr+th+grouping.csv"
#dir1="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/identification/"
#dir2="/data/SBCS-BessantLab/shyama/Data/Bristol/Signature_peptide_seeded_database/Mosquito/"
#host="Aedes aegypti"
#virus="_MCFA"

Mats=readMat(f1,f2,dir1,dir2,rev,peptide,pepThreshold,upper)
peptides=readPSM(f3,dir1)
sink(file="/data/home/btw796/Code2/Proteomics/R/proViralLogBatHendraReviewed.txt",append=FALSE)
#host="HUMAN"
#virus="_ADE05"
trinityORFStart="^comp"
print("Total PAG count:")
print(length(unique(Mats[[1]][,'PAG.ID'])))
print("Known Virus Count(Unique|Total):")
print(length(unique(Mats[[1]][grep(virus,Mats[[1]][,'protein.accession']),'protein.accession'])))
print(length(Mats[[1]][grep(virus,Mats[[1]][,'protein.accession']),'protein.accession']))

print("Known Virus Count PAG(Unque|Total)")
print(length(unique(Mats[[1]][grep(virus,Mats[[1]][,'protein.accession']),'PAG.ID'])))
print(length(Mats[[1]][grep(virus,Mats[[1]][,'protein.accession']),'PAG.ID']))
vProt=unique(Mats[[1]][grep(virus,Mats[[1]][,'protein.accession']),'description'])
print("Known Virus Peptide Count:")
print(length(unique(peptideCount(peptides,vProt))))
vPAG=unique(Mats[[1]][grep(virus,Mats[[1]][,'protein.accession']),'PAG.ID'])

print("Host Count (Unique|Total):")
print(length(unique(Mats[[1]][grep(host,Mats[[1]][,'protein.accession']),'protein.accession'])))
print(length(Mats[[1]][grep(host,Mats[[1]][,'protein.accession']),'protein.accession']))

print("Host Count PAG(Unique|Total):")
print(length(unique(Mats[[1]][grep(host,Mats[[1]][,'protein.accession']),'PAG.ID'])))
print(length(Mats[[1]][grep(host,Mats[[1]][,'protein.accession']),'PAG.ID']))
hMat=as.data.frame(Mats[[1]][grep(host,Mats[[1]][,'protein.accession']),c("protein.accession","PAG.ID")])
print(dim(hMat))
hMat=hMat[!duplicated(hMat[,c("protein.accession","PAG.ID")]),]
print(dim(hMat))
PAGCount(hMat)
hProt=unique(Mats[[1]][grep(host,Mats[[1]][,'protein.accession']),'description'])
print("Known Host Peptide Count:")
print(length(unique(peptideCount(peptides,hProt))))

hPAG=unique(Mats[[1]][grep(host,Mats[[1]][,'protein.accession']),'PAG.ID'])

print("Unknown Virus Count(Unique|Total):")
print(length(unique(Mats[[1]][-grep(paste(host,virus,trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'protein.accession'])))
print(length(Mats[[1]][-grep(paste(host,virus,trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'protein.accession']))
print(unique(Mats[[1]][-grep(paste(host,virus,trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'protein.accession']))

print("Unknown Virus Count PAG(Unique|Total):")
print(length(unique(Mats[[1]][-grep(paste(host,virus,trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'PAG.ID'])))
print(length(Mats[[1]][-grep(paste(host,virus,trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'PAG.ID']))
uPAG=unique(Mats[[1]][-grep(paste(host,virus,trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'PAG.ID'])
print(uPAG)
uProt=unique(Mats[[1]][-grep(paste(host,virus,trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'description'])
print("Known Unknown Virus Peptide Count:")
print(length(unique(peptideCount(peptides,uProt))))
print("No Match:")
print(length(unique(Mats[[1]][grep(paste(trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'protein.accession'])))
print("No Match PAG:")
print(length(unique(Mats[[1]][grep(paste(trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'PAG.ID'])))

nPAG=unique(Mats[[1]][grep(paste(trinityORFStart,sep="|"),Mats[[1]][,'protein.accession']),'PAG.ID'])

print("Overlaping PAG Host vs Virus")
print(length(which(hPAG %in% vPAG)))

print("Overlaping PAG Host vs Unknown Virus")
print(length(which(hPAG %in% uPAG)))

print("Overlaping PAG Virus vs Unknown Virus")
print(length(which(vPAG %in% uPAG)))

print("Overlaping PAG Host vs no Match")
print(length(which(hPAG %in% nPAG)))

print("Overlaping PAG Virus vs no Match")
print(length(which(vPAG %in% nPAG)))

print("Overlaping PAG unknown virus vs no Match")
print(length(which(uPAG %in% nPAG)))

print("Identified ORFs Count:")
print(dim(Mats[[1]]))
print("BLAST DB:")
print(dim(Mats[[2]]))
print("Unique DB:")
print(length(unique(Mats[[2]][,'hit_def'])))
print("Host DB(Unique,ORFs count):")
print(length(unique(Mats[[2]][grep(host,Mats[[2]][,'hit_def']),'hit_def'])))
print(length(Mats[[2]][grep(host,Mats[[2]][,'hit_def']),'hit_def']))

print("Known Virus DB(Unique,ORFs count):")
print(length(unique(Mats[[2]][grep(virus,Mats[[2]][,'hit_def']),'hit_def'])))
print(length(Mats[[2]][grep(virus,Mats[[2]][,'hit_def']),'hit_def']))

print("Unknown Virus DB(Unique,ORFs count):")
print(length(unique(Mats[[2]][-grep(paste(host,virus,sep="|"),Mats[[2]][,'hit_def']),'hit_def'])))
print(length(Mats[[2]][-grep(paste(host,virus,sep="|"),Mats[[2]][,'hit_def']),'hit_def']))
sink()

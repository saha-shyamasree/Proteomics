readExclusiveProteinList <- function(filename, dir)
{
    as.matrix(read.table(file=paste(dir,filename,sep=""), sep='\t',header=TRUE,quote = ""))
}

readBioMartMapping <- function(filename, dir)
{
    as.matrix(read.csv(file=paste(dir,filename,sep=""),header=TRUE, quote = ""))
}

uniIdRegex<-function(Mat,column)
{
    Mat[,column]=sub("(.*)?\\|([^\\|]+)\\|(.*)?","\\2",Mat[,column])
    Mat
}

blast<-function(filename,folder)
{
    read.csv(file=paste(folder,filename,sep=""),header=TRUE,quote = "")
}

blastFilterMatch <- function(blastMat,flag)
{
    if(flag==1)
    {
        blastMat[which(blastMat[,'match']=='yes'),]
    }
    else
    {
        if(flag==0)
        {
            blastMat[which(blastMat[,'match']=='no'),]
        }
        else
        {
            NULL;
        }
    }
}

blastFilterEval <- function(blastMat,upper)
{
    print(paste("evalue=",upper))
    subset(blastMat,e.value<=upper)
}

isoformCheck <- function(id, Mat)
{
    
}

excUni="onlyDM_from_raw_uniprot+fdr+th+grouping+prt.csv1e-30_compared_trinity_PITORF+fdr+th+grouping+prt.csv_1e-30.tsv"

##################################################################### Finding Uniprot protein that can not be identified from trinity ORFs, because they were not expressed  #################################
excUni="PAGComparison_UniprotOnly_withTrinity.tsv"
d1="D:/data/Results/Human-Adeno/Identification/"

uniOnly=readExclusiveProteinList(excUni, d1)
uniOnly_ids=uniIdRegex(uniOnly,3)

biomart="GRCh38_gene_transcript_protein_biomart.csv"
d4="D:/data/Data/"
biomartMat=readBioMartMapping(biomart,d4)
d3="D:/data/blast/blastCSV/"
f4="human_adenovirus_trinityRNA_blast.csv"
upper=0.000000000000000000000000000001
blastDB=blastFilterEval(blastFilterMatch(blast(f4,d3),1),upper)
f5="AnchortoAnchor_PAGs_uni_trinity1e-30.tsv"
d2="D:/data/Results/Human-Adeno/Identification/"
overlapMat=readExclusiveProteinList(f5,d2)
f6="AnchortoSub_PAGs_uni_trinity1e-30.tsv"
overlapMat=rbind(overlapMat,readExclusiveProteinList(f6,d2))
uni_ids=sub("([^-]+)?-(.*)?","\\1",uniOnly_ids[,'protein.accession'])
overlap_uni_ids=sub("(.*)?\\|([^\\|]+)?\\|(.*)?","\\2",overlapMat[,'protein.accession'])
length(which(uni_ids %in% overlap_uni_ids))
iso_overlap=uni_ids[which(uni_ids %in% overlap_uni_ids)]

###PAG level overlap
### Finding ensembl transcript ids: no of transcripts can be more as one protein might be produced from different transcript
uni_ids_backup=uni_ids
uni_ids=uni_ids[-which(uni_ids %in% overlap_uni_ids)]
length(uni_ids)

transcripts_ids=biomartMat[which(biomartMat[,'UniProt.SwissProt.Accession'] %in% uni_ids),]
rest_uni_ids=uni_ids[-which(uni_ids %in% biomartMat[,'UniProt.SwissProt.Accession'])]
biomartMat_temp=biomartMat[which(!biomartMat[,'UniProt.SwissProt.Accession'] %in% uni_ids),]
trembl_transcripts_ids=biomartMat_temp[which(biomartMat_temp[,'UniProt.TrEMBL.Accession'] %in% rest_uni_ids),]
blastTranscripts=sub("([^\\s]+) (.*)?","\\1",blastDB[,'hit_def'])
rest2_uni_ids=rest_uni_ids[-which(rest_uni_ids %in% biomartMat_temp[,'UniProt.TrEMBL.Accession'])]
#transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.SwissProt.Accession']

UniprotOnlyBioMart=rbind(transcripts_ids,trembl_transcripts_ids)

length(unique(transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.SwissProt.Accession']))
length(unique(trembl_transcripts_ids[which(trembl_transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.TrEMBL.Accession']))
expressedUniProteins=c(unique(transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.SwissProt.Accession']),unique(trembl_transcripts_ids[which(trembl_transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.TrEMBL.Accession']))
old_proteins_trinity=uni_ids[which(!uni_ids %in% expressedUniProteins)]

d5="D:/data/Results/Human-Adeno/GIOPaperResults/PostProcessings/"
write.csv(old_protein,file=paste(d5,"Unexpressed_trinity.csv",sep=""),quote=FALSE,row.names = FALSE)
################################################################################################################

##################################################################### Finding Uniprot protein that can not be identified from cufflinks ORFs, because they were not expressed  #################################
excUni="PAGComparison_UniprotOnly_withCuff1e-30.tsv"
d1="D:/data/Results/Human-Adeno/Identification/"

uniOnly=readExclusiveProteinList(excUni, d1)
uniOnly_ids=uniIdRegex(uniOnly,3)

biomart="GRCh38_gene_transcript_protein_biomart.csv"
d4="D:/data/Data/"
biomartMat=readBioMartMapping(biomart,d4)
d3="D:/data/blast/blastCSV/"
f4="human_CufflinksRNA_blast.csv"
upper=0.000000000000000000000000000001
blastDB=blastFilterEval(blastFilterMatch(blast(f4,d3),1),upper)
f5="AnchortoAnchor_PAGs_uni_cuff1e-30.tsv"
d2="D:/data/Results/Human-Adeno/Identification/"
overlapMat=readExclusiveProteinList(f5,d2)
f6="AnchortoSub_PAGs_uni_cuff1e-30.tsv"
overlapMat=rbind(overlapMat,readExclusiveProteinList(f6,d2))
uni_ids=sub("([^-]+)?-(.*)?","\\1",uniOnly_ids[,'protein.accession'])
overlap_uni_ids=sub("(.*)?\\|([^\\|]+)?\\|(.*)?","\\2",overlapMat[,'protein.accession'])
length(which(uni_ids %in% overlap_uni_ids))
iso_overlap=uni_ids[which(uni_ids %in% overlap_uni_ids)]

###PAG level overlap
### Finding ensembl transcript ids: no of transcripts can be more as one protein might be produced from different transcript
uni_ids_backup=uni_ids
uni_ids=uni_ids[-which(uni_ids %in% overlap_uni_ids)]
length(uni_ids)

transcripts_ids=biomartMat[which(biomartMat[,'UniProt.SwissProt.Accession'] %in% uni_ids),]
rest_uni_ids=uni_ids[-which(uni_ids %in% biomartMat[,'UniProt.SwissProt.Accession'])]
biomartMat_temp=biomartMat[which(!biomartMat[,'UniProt.SwissProt.Accession'] %in% uni_ids),]
trembl_transcripts_ids=biomartMat_temp[which(biomartMat_temp[,'UniProt.TrEMBL.Accession'] %in% rest_uni_ids),]
blastTranscripts=sub("([^\\s]+) (.*)?","\\1",blastDB[,'hit_def'])
rest2_uni_ids=rest_uni_ids[-which(rest_uni_ids %in% biomartMat_temp[,'UniProt.TrEMBL.Accession'])]
#transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.SwissProt.Accession']

UniprotOnlyBioMart=rbind(transcripts_ids,trembl_transcripts_ids)

length(unique(transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.SwissProt.Accession']))
length(unique(trembl_transcripts_ids[which(trembl_transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.TrEMBL.Accession']))
expressedUniProteins=c(unique(transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.SwissProt.Accession']),unique(trembl_transcripts_ids[which(trembl_transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'UniProt.TrEMBL.Accession']))
old_proteins_cufflinks=uni_ids[which(!uni_ids %in% expressedUniProteins)]

d5="D:/data/Results/Human-Adeno/GIOPaperResults/PostProcessings/"
write.csv(old_protein,file=paste(d5,"Unexpressed_cufflinks.csv",sep=""),quote=FALSE,row.names = FALSE)
################################################################################################################
## Overlap between trinity and cufflinks old protein
################################################################################################################
length(which(old_proteins_trinity %in% old_proteins_cufflinks))



################################################################################################################

#unique(transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),'Ensembl.Transcript.ID'])
##Ensembl transcripts for Uni only proteins, which have Ensemble transcripts to uniprot/trembl protein map and Trinity transcript mapped to them.

uniprot_corresponding_transcripts=unique(transcripts_ids[which(transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),])
uniprot_corresponding_trembl_transcripts=unique(trembl_transcripts_ids[which(trembl_transcripts_ids[,'Ensembl.Transcript.ID'] %in% unique(blastTranscripts)),])
uniprot_corresponding_transcripts_all=rbind(uniprot_corresponding_transcripts,uniprot_corresponding_trembl_transcripts)

trinityTranscripts_all=lapply(uniprot_corresponding_transcripts_all[,2],function(x,y){grep(paste("^",x,sep=""),y[,"hit_def"])},y=blastDB)

trinityTranscriptsIDX_all=unlist(trinityTranscripts_all)
##subset of blastDB
UniTranscriptfromTrinity_all = blastDB[trinityTranscriptsIDX_all,]
blastDB=cbind(blastDB,sub("(.*)?FPKM=([^\\s]+)?\\s(.*)?","\\2",blastDB[,'query_name']))
fpkm=sub("(.*)?FPKM=([^\\s]+)?\\s(.*)?","\\2",UniTranscriptfromTrinity_all[,'query_name'])
UniTranscriptfromTrinity_all=cbind(UniTranscriptfromTrinity_all,fpkm)
UniTranscriptfromTrinity_all[,22]=as.numeric(as.character(UniTranscriptfromTrinity_all[,22]))

UniTranscriptfromTrinity_all=cbind(UniTranscriptfromTrinity_all,(UniTranscriptfromTrinity_all[,'query_length']/UniTranscriptfromTrinity_all[,'hit_length']))

UniTranscriptFPKM50Plus=UniTranscriptfromTrinity_all[which(UniTranscriptfromTrinity_all[,22]>50),]
UniTranscriptFPKMbelowOne=UniTranscriptfromTrinity_all[which(UniTranscriptfromTrinity_all[,22]<1),]

unique(as.character(UniTranscriptfromTrinity_all[which(UniTranscriptfromTrinity_all[,'long_match']>=0.5),'hit_def']))
dim(unique(UniTranscriptfromTrinity_all))
UniTranscriptfromTrinity_all_unique=unique(UniTranscriptfromTrinity_all)

transcriptLevel=sub("([^\\s]+?)\\s(.*)?","\\1",UniTranscriptfromTrinity_all_unique[,'query_name'])
indxs=unlist(lapply(transcriptLevel,function(x,y){grep(x,y[,'query_name'])},y=blastDBPrt))
all_trans_fpkm=sub("(.*)?FPKM=([^\\s]+)?\\s(.*)?","\\2",blastDB[,'query_name'])
all_trans_fpkm100Below= all_trans_fpkm[which( all_trans_fpkm<100)]

##################### Overlapping proteins #######################

overlapAtoAIDS=sub("(.*)?\\|([^\\|]+)?\\|(.*)?","\\2",a_to_a[,'protein.accession'])
trinityAnchorOverlapIDS=sub("(.*)?\\|([^\\|]+)?\\|(.*)?","\\2",m1_without_atoaPAG_anchor[res$idx,'protein.accession'])
uniAnchorOverlapIDS=sub("(.*)?\\|([^\\|]+)?\\|(.*)?","\\2",Mat2_without_subtoa_anchor[res2$idx,'protein.accession'])
overlapAnchorIDS=c(overlapAtoAIDS,trinityAnchorOverlapIDS,uniAnchorOverlapIDS)

tids=which(biomartMat[,'UniProt.SwissProt.Accession'] %in% overlapAnchorIDS)
overlap_transcripts=biomartMat[tids,]
biomart_overlap_temp=biomartMat[-tids,]
overlap_transcripts=rbind(overlap_transcripts,biomart_overlap_temp[which(biomart_overlap_temp[,'UniProt.TrEMBL.Accession'] %in% overlapAnchorIDS),])
unique_transcripts=unique(overlap_transcripts[,'Ensembl.Transcript.ID'])

overlaptrinityTranscripts_all=lapply(unique_transcripts,function(x,y){grep(paste("^",x,sep=""),y[,"hit_def"])},y=blastDB)
overlaptrinityTranscriptsIDX_all=unlist(overlaptrinityTranscripts_all)
overlapTranscript= blastDB[overlaptrinityTranscriptsIDX_all,]

fpkm=sub("(.*)?FPKM=([^\\s]+)?\\s(.*)?","\\2",overlapTranscript[,'query_name'])
overlapTranscript=cbind(overlapTranscript,fpkm)

overlapTranscript[,'fpkm']=as.numeric(overlapTranscript[,'fpkm'])



f6="trinity_PITORF_human_adeno_blast2.csv"
upper=0.000000000000000000000000000001
blastDBProteins=blastFilterEval(blastFilterMatch(blast(f6,d3),1),upper)




exp_ids=which(blastTranscripts %in% unique(transcripts_ids[,'Ensembl.Transcript.ID']))
exp_ids2=which(blastTranscripts %in% unique(trembl_transcripts_ids[,'Ensembl.Transcript.ID']))
exp_ids_all=c(exp_ids,exp_ids2)

uni_transcripts = blastDB[exp_ids_all,]





#overlap_ids=sub("(.*)?\\|([^\\|]+)?\\|(.*)?","\\2",overlapMat[,'protein.accession'])

overlap_ids_no_iso=sub("([^-]+)?-(.*)?","\\1",overlap_ids)
uniq_overlap_no_iso=unique(overlap_ids_no_iso)
uniq_uni_ids=unique(uni_ids)
uni_ids_for_biomart_check=uni_ids[which(!uni_ids %in% uniq_overlap_no_iso)]
length(which(uni_ids_for_biomart_check %in% biomartMat[,'UniProt.SwissProt.Accession']))
[1] 236
length(which(uni_ids_for_biomart_check %in% biomartMat[,'UniProt.TrEMBL.Accession']))
[1] 279
length(unique(biomartMat[which(biomartMat[,'UniProt.SwissProt.Accession'] %in% uni_ids_for_biomart_check),'Ensembl.Transcript.ID']))

gene_transcripts_ids=biomartMat[which(biomartMat[,'UniProt.SwissProt.Accession'] %in% uni_ids_for_biomart_check),c('Ensembl.Gene.ID','Ensembl.Transcript.ID','UniProt.SwissProt.Accession')]
biomartMat_temp=biomartMat[which(!biomartMat[,'UniProt.SwissProt.Accession'] %in% uni_ids_for_biomart_check),]
trembl_gene_transcripts_ids=biomartMat_temp[which(biomartMat_temp[,'UniProt.TrEMBL.Accession'] %in% uni_ids_for_biomart_check),c('Ensembl.Gene.ID','Ensembl.Transcript.ID','UniProt.TrEMBL.Accession')]

blastTranscripts=sub("([^\\s]+) (.*)?","\\1",blastDB[,'hit_def'])
blastGene=sub("(.*)?gene:([^\\s]+) (.*)?","\\2",blastDB[,'hit_def'])

length(which(unique(blastTranscripts) %in% unique(gene_transcripts_ids[,'Ensembl.Transcript.ID'])))


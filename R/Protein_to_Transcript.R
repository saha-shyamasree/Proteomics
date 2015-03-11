
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

getTranscriptTypes<-function(Mat)
{
    transcriptBioTypes=sub("(.*) transcript_biotype:(.*)?","\\2",Mat[,'hit_def'])
    transcriptBioTypes
}

addTranscriptTypeColumn <- function(Mat)
{
    Mat=cbind(Mat,getTranscriptTypes(Mat))
    Mat
}

blastFilterEval <- function(blastMat,upper)
{
    print(paste("evalue=",upper))
    subset(blastMat,e.value<=upper)
}

blastFilterTranscriptType<-function(Mat,type,column)
{
    Mat[which(Mat[,column]==type),]
}

queryIDs <- function(Mat)
{
    sub("([^\\s]+)? (.*)?","\\1",Mat[,'query_name'])
}

proteinIDs <- function(Mat)
{
    strsplit(Mat,",",fixed=TRUE)
}

removeORFid <- function(proteinIds)
{
    sub("(.*)?_ORF(.*)?","\\1",proteinIds)
}

removeuORFid <- function(proteinIds)
{
    sub("(.*)?_u(.*)?","\\1",proteinIds)
}

removeuORFidFromList <- function(proteinList)
{
    subsList=list()
    for(i in 1:length(proteinList))
    {
        subsList[[i]]=sub("(.*)?_u(.*)?","\\1",proteinList[[i]])
    }
    subsList
}

removeORFidFromList <- function(proteinList)
{
    subsList=list()
    for(i in 1:length(proteinList))
    {
        subsList[[i]]=sub("(.*)?_ORF(.*)?","\\1",proteinList[[i]])
    }
    subsList
}



groupProteins<- function(Mat)
{
    DBA=Mat[grep("^Dataset_A_",Mat[,'protein.accession']),]
    DBB=Mat[grep("^Dataset_B_",Mat[,'protein.accession']),]
    DBC=Mat[grep("^Dataset_C_",Mat[,'protein.accession']),]
    list(A=DBA,B=DBB,C=DBC)
}

removeDBprefix<-function(Mat)
{
    sub("(.*)?_(comp(.*)?)","\\2",Mat)
}
proteinGroup <- function(filename,dirc)
{
    as.matrix(read.csv(file=paste(dirc,filename,sep=""),header=TRUE))
}

proteinGroupFiltered<- function(proteinGrp,rev,peptide,pepThreshold)
{
    proteinGrp=as.matrix(proteinGrp)
    print(dim(proteinGrp))
    if(rev==1)
    {
        proteinGrp=proteinGrp[-(grep("_REVERSED",proteinGrp[,'protein.accession'])),]
        print(dim(proteinGrp))
    }
    if(peptide==1)
    {
        proteinGrp[,'distinct.peptide.sequences']=as.numeric(proteinGrp[,'distinct.peptide.sequences'])
        proteinGrp=proteinGrp[which(proteinGrp[,'distinct.peptide.sequences']>pepThreshold),]
        print(dim(proteinGrp))
    }
    proteinGrp
}

trV5Pep="trinity_PITORF+fdr+th+groupingV5.csv"
trV5="trinity_PITORF+fdr+th+grouping+prtV5.csv"
trD="C:/Users/shyama/Dropbox/results/Human_adenovirus/sORF/"
rnaBlast="human_adenovirus_trinityRNA_blast.csv"
blastDBDir="D:/data/blast/blastCSV/"
blastDB=addTranscriptTypeColumn(blastFilterEval(blastFilterMatch(blast(rnaBlast,blastDBDir),1),upper)) #//blast(rnaBlast,blastDBDir)

Mat1=proteinGroupFiltered(proteinGroup(trV5,trD),rev,peptide,pepThreshold)
Mat1Pep=proteinGroup(trV5Pep,trD)


ORFGrp=groupProteins(Mat1)
orfIdC=removeuORFid(unlist(proteinIDs(removeDBprefix(ORFGrp$C[,'protein.accession']))))
orfIdA=removeORFid(unlist(proteinIDs(removeDBprefix(ORFGrp$A[,'protein.accession']))))
orfIdB=removeuORFid(unlist(proteinIDs(removeDBprefix(ORFGrp$B[,'protein.accession']))))

length(ORFGrp$C[,'protein.accession'])
length(ORFGrp$B[,'protein.accession'])
length(ORFGrp$A[,'protein.accession'])



DBA=removeORFidFromList(proteinIDs(removeDBprefix(ORFGrp$A[,'protein.accession'])))
DBB=removeuORFidFromList(proteinIDs(removeDBprefix(ORFGrp$B[,'protein.accession'])))
DBC=removeuORFidFromList(proteinIDs(removeDBprefix(ORFGrp$C[,'protein.accession'])))

blasDBCDS=blastFilterTranscriptType(blastDB,'protein_coding',22)
CDSIDs=queryIDs(blasDBCDS)
length(which(CDSIDs %in% orfIdA))
length(which(CDSIDs %in% orfIdB))
length(which(CDSIDs %in% orfIdC))

DBA_CDSIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=CDSIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_CDSIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_CDSIDs_length,z=DBA_length)))
length(unlist(lapply(which(DBA_CDSIDs_length>0),function(x,y,z){if(y[x]!=z[x]){x}},y=DBA_CDSIDs_length,z=DBA_length)))

blasDBnonsense_mediated_decay=blastFilterTranscriptType(blastDB,'nonsense_mediated_decay',22)
nonsense_mediated_decayIDs=queryIDs(blasDBnonsense_mediated_decay)
length(nonsense_mediated_decayIDs)
length(which(nonsense_mediated_decayIDs %in% orfIdA))
length(which(nonsense_mediated_decayIDs %in% orfIdB))
length(which(nonsense_mediated_decayIDs %in% orfIdC))

DBA_nonsense_mediated_decayIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=nonsense_mediated_decayIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_nonsense_mediated_decayIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_nonsense_mediated_decayIDs_length,z=DBA_length)))
length(unlist(lapply(which(DBA_nonsense_mediated_decayIDs_length>0),function(x,y,z){if(y[x]!=z[x]){x}},y=DBA_nonsense_mediated_decayIDs_length,z=DBA_length)))

blasDBLinc=blastFilterTranscriptType(blastDB,'lincRNA',22)
lincIDs=queryIDs(blasDBLinc)
length(lincIDs)
length(which(lincIDs %in% orfIdA))
length(which(lincIDs %in% orfIdB))
length(which(lincIDs %in% orfIdC))

DBA_lincIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=lincIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_lincIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_lincIDs_length,z=DBA_length)))
length(unlist(lapply(which(DBA_lincIDs_length>0),function(x,y,z){if(y[x]!=z[x]){x}},y=DBA_lincIDs_length,z=DBA_length)))

grporfnIds=DBA[unlist(lapply(which(DBA_nonsense_mediated_decayIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_nonsense_mediated_decayIDs_length,z=DBA_length))]
idx=c()
idxPep=c()
for(i in 1:length(grporfnIds))
{
    j=grep(paste("^Dataset_A_",grporfnIds[[i]][1],sep=""),Mat1[,'protein.accession'])
    idx=c(idx,j)
    idxPep=c(idxPep,grep(paste("(.*)?",Mat1[j,'protein.accession'],"(.*)?",sep=""), Mat1Pep[,'proteinacc_start_stop_pre_post_.']))
}
Mat1_nonsense_mediated_decay=Mat1[idx,]

###following way of counting peptide is not right because many of these ORFs are sub-member of the PAG, hence peptides associated to it is not found from the row.
ri_unq_pep=strsplit(Mat1_nonsense_mediated_decay[which(Mat1_nonsense_mediated_decay[,'unique.peptides']!=""),'unique.peptides'],';')

ri_raz_pep=strsplit(Mat1_nonsense_mediated_decay[which(Mat1_reatainedIntron[,'razor.peptides']!=""),'razor.peptides'],';')
ri_peps=sub("^\\s+", "", unlist(c(ri_unq_pep,ri_raz_pep)))
Mat1Pepnonsense_mediated_decayPep=Mat1Pepnonsense_mediated_decay[which(Mat1Pepnonsense_mediated_decay[,'Sequence'] %in% ri_peps),]

#####################################

Mat1Pepnonsense_mediated_decay=Mat1Pep[idxPep,]
dim(Mat1Pep[which(Mat1Pep[,'Sequence'] %in% ri_peps),])
rorfIds=unlist(DBA[unlist(lapply(which(DBA_nonsense_mediated_decayIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_nonsense_mediated_decayIDs_length,z=DBA_length))])##all of the transcript ids from the ORF id matches to a 'processed_transcript'.
which(nonsense_mediated_decayIDs %in% rorfIds)
write.table(blasDBnonsense_mediated_decay[which(nonsense_mediated_decayIDs %in% rorfIds),],file=paste(blastDBDir,"nonsense_mediated_decay.csv",sep=""),sep=",",quote = FALSE,row.names = FALSE)
write.table(Mat1[idx,],file=paste(trD,"nonsense_mediated_decayPrt.tsv",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
write.table(Mat1Pep[idxPep,],file=paste(trD,"nonsense_mediated_decay.tsv",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
################################################################################################################################


blasDBretained_intron=blastFilterTranscriptType(blastDB,'retained_intron',22)
retained_intronIDs=queryIDs(blasDBretained_intron)
length(retained_intronIDs)
length(which(retained_intronIDs %in% orfIdA))
length(which(retained_intronIDs %in% orfIdB))
length(which(retained_intronIDs %in% orfIdC))

DBA_retained_intronIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=retained_intronIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_retained_intronIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_retained_intronIDs_length,z=DBA_length)))
length(unlist(lapply(which(DBA_retained_intronIDs_length>0),function(x,y,z){if(y[x]!=z[x]){x}},y=DBA_retained_intronIDs_length,z=DBA_length)))

grprorfIds=DBA[unlist(lapply(which(DBA_retained_intronIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_retained_intronIDs_length,z=DBA_length))]
idx1=c()
idxPep1=c()
for(i in 1:length(grprorfIds))
{
    j=grep(paste("^Dataset_A_",grprorfIds[[i]][1],sep=""),Mat1[,'protein.accession'])
    idx1=c(idx1,j)
    idxPep1=c(idxPep1,grep(paste("(.*)?",Mat1[j,'protein.accession'],"(.*)?",sep=""), Mat1Pep[,'proteinacc_start_stop_pre_post_.']))
}
Mat1_reatainedIntron=Mat1[idx,]

###following way of counting peptide is not right because many of these ORFs are sub-member of the PAG, hence peptides associated to it is not found from the row.
ri_unq_pep=strsplit(Mat1_reatainedIntron[which(Mat1_reatainedIntron[,'unique.peptides']!=""),'unique.peptides'],';')

ri_raz_pep=strsplit(Mat1_reatainedIntron[which(Mat1_reatainedIntron[,'razor.peptides']!=""),'razor.peptides'],';')
ri_peps=sub("^\\s+", "", unlist(c(ri_unq_pep,ri_raz_pep)))
#####################################

Mat1PepRetainedIntron=Mat1Pep[idxPep,]
Mat1PepRetainedIntronPep=Mat1PepRetainedIntron[which(Mat1PepRetainedIntron[,'Sequence'] %in% ri_peps),]
dim(Mat1Pep[which(Mat1Pep[,'Sequence'] %in% ri_peps),])
rorfIds=unlist(DBA[unlist(lapply(which(DBA_retained_intronIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_retained_intronIDs_length,z=DBA_length))])##all of the transcript ids from the ORF id matches to a 'processed_transcript'.
which(retained_intronIDs %in% rorfIds)
write.table(blasDBretained_intron[which(retained_intronIDs %in% rorfIds),],file=paste(blastDBDir,"retainedIntron.csv",sep=""),sep=",",quote = FALSE,row.names = FALSE)
write.table(Mat1[idx,],file=paste(trD,"retained_intronPrt.tsv",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
write.table(Mat1Pep[idxPep,],file=paste(trD,"retained_intron.tsv",sep=""),sep="\t",quote = FALSE,row.names = FALSE)

which(DBA_retained_intronIDs_length>0) ###at least one of the transcript ids from the ORF id matches to a 'processed_transcript'.

#Dataset_A_comp2526_c1_seq3_ORF9_Frame_3_300-965_102_112_R_T
#Dataset_A_comp2526_c1_seq3_ORF9_Frame_3_300-965_159_171_R_R

blasDBantisense=blastFilterTranscriptType(blastDB,'antisense',22)
antisenseIDs=queryIDs(blasDBantisense)
length(antisenseIDs)
length(which(antisenseIDs %in% orfIdA))
length(which(antisenseIDs %in% orfIdB))
length(which(antisenseIDs %in% orfIdC))

DBA_antisenseIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=antisenseIDs),length))
DBA_length=unlist(lapply(DBA,length))
##all parents are from the same bio-type
length(unlist(lapply(which(DBA_antisenseIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_antisenseIDs_length,z=DBA_length)))
##not-all parents are from the same bio-type
length(unlist(lapply(which(DBA_antisenseIDs_length>0),function(x,y,z){if(y[x]!=z[x]){x}},y=DBA_antisenseIDs_length,z=DBA_length)))


blasDBprocessed_transcript=blastFilterTranscriptType(blastDB,'processed_transcript',22)
processed_transcriptIDs=queryIDs(blasDBprocessed_transcript)
length(processed_transcriptIDs)
length(which(processed_transcriptIDs %in% orfIdA))
length(which(processed_transcriptIDs %in% orfIdB))
length(which(processed_transcriptIDs %in% orfIdC))

DBA_processed_transcriptIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=processed_transcriptIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_processed_transcriptIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_processed_transcriptIDs_length,z=DBA_length)))
length(unlist(lapply(which(DBA_processed_transcriptIDs_length>0),function(x,y,z){if(y[x]!=z[x]){x}},y=DBA_processed_transcriptIDs_length,z=DBA_length)))

processed_transcriptIdentified=processed_transcriptIDs[which(processed_transcriptIDs %in% orfIdA)]
processed_transcriptIdentified=c(processed_transcriptIdentified,processed_transcriptIDs[which(processed_transcriptIDs %in% orfIdB)])
processed_transcriptIdentified=c(processed_transcriptIdentified,processed_transcriptIDs[which(processed_transcriptIDs %in% orfIdC)])
length(processed_transcriptIdentified)
length(unique(processed_transcriptIdentified))
unq_processed_transcriptIdentified=unique(processed_transcriptIdentified)


DBA[unlist(lapply(which(DBA_processed_transcriptIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_processed_transcriptIDs_length,z=DBA_length))]##all of the transcript ids from the ORF id matches to a 'processed_transcript'.

which(DBA_processed_transcriptIDs_length>0) ###at least one of the transcript ids from the ORF id matches to a 'processed_transcript'.

DBB_processed_transcriptIDs_length = unlist(lapply(lapply(DBB,function(x,y){which(x %in% y)}, y=processed_transcriptIDs),length))
DBB_length=unlist(lapply(DBB,length))

DBB[unlist(lapply(which(DBB_processed_transcriptIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBB_processed_transcriptIDs_length,z=DBB_length))] ##all of the transcript ids from the ORF id matches to a 'processed_transcript'.
which(DBB_processed_transcriptIDs_length>0) ###at least one of the transcript ids from the ORF id matches to a 'processed_transcript'.

blasDBMt_rRNA=blastFilterTranscriptType(blastDB,'Mt_rRNA',22)
Mt_rRNAIDs=queryIDs(blasDBMt_rRNA)
length(Mt_rRNAIDs)
length(which(Mt_rRNAIDs %in% orfIdA))
length(which(Mt_rRNAIDs %in% orfIdB))
length(which(Mt_rRNAIDs %in% orfIdC))

DBA_Mt_rRNAIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=Mt_rRNAIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_Mt_rRNAIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_Mt_rRNAIDs_length,z=DBA_length)))
length(unlist(lapply(which(DBA_Mt_rRNAIDs_length>0),function(x,y,z){if(y[x]!=z[x]){x}},y=DBA_Mt_rRNAIDs_length,z=DBA_length)))

blasDBprocessed_pseudogene=blastFilterTranscriptType(blastDB,'processed_pseudogene',22)
processed_pseudogeneIDs=queryIDs(blasDBprocessed_pseudogene)
length(processed_pseudogeneIDs)
length(which(processed_pseudogeneIDs %in% orfIdA))
length(which(processed_pseudogeneIDs %in% orfIdB))
length(which(processed_pseudogeneIDs %in% orfIdC))

DBA_processed_pseudogeneIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=processed_pseudogeneIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_processed_pseudogeneIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_processed_pseudogeneIDs_length,z=DBA_length)))


blasDBmisc_RNA=blastFilterTranscriptType(blastDB,'misc_RNA',22)
misc_RNAIDs=queryIDs(blasDBmisc_RNA)
length(misc_RNAIDs)
length(which(misc_RNAIDs %in% orfIdA))
length(which(misc_RNAIDs %in% orfIdB))
length(which(misc_RNAIDs %in% orfIdC))

DBA_misc_RNAIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=misc_RNAIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_misc_RNAIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_misc_RNAIDs_length,z=DBA_length)))


blasDBMt_tRNA=blastFilterTranscriptType(blastDB,'Mt_tRNA',22)
Mt_tRNAIDs=queryIDs(blasDBMt_tRNA)
length(Mt_tRNAIDs)
length(which(Mt_tRNAIDs %in% orfIdA))
length(which(Mt_tRNAIDs %in% orfIdB))
length(which(Mt_tRNAIDs %in% orfIdC))

DBA_Mt_tRNAIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=Mt_tRNAIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_Mt_tRNAIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_Mt_tRNAIDs_length,z=DBA_length)))


blasDBtranscribed_unprocessed_pseudogene=blastFilterTranscriptType(blastDB,'transcribed_unprocessed_pseudogene',22)
transcribed_unprocessed_pseudogeneIDs=queryIDs(blasDBtranscribed_unprocessed_pseudogene)
length(transcribed_unprocessed_pseudogeneIDs)
length(which(transcribed_unprocessed_pseudogeneIDs %in% orfIdA))
length(which(transcribed_unprocessed_pseudogeneIDs %in% orfIdB))
length(which(transcribed_unprocessed_pseudogeneIDs %in% orfIdC))

DBA_transcribed_unprocessed_pseudogeneIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=transcribed_unprocessed_pseudogeneIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_transcribed_unprocessed_pseudogeneIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_transcribed_unprocessed_pseudogeneIDs_length,z=DBA_length)))


blasDBunprocessed_pseudogene=blastFilterTranscriptType(blastDB,'unprocessed_pseudogene',22)
unprocessed_pseudogeneIDs=queryIDs(blasDBunprocessed_pseudogene)
length(unprocessed_pseudogeneIDs)
length(which(unprocessed_pseudogeneIDs %in% orfIdA))
length(which(unprocessed_pseudogeneIDs %in% orfIdB))
length(which(unprocessed_pseudogeneIDs %in% orfIdC))

DBA_unprocessed_pseudogeneIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=unprocessed_pseudogeneIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_unprocessed_pseudogeneIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_unprocessed_pseudogeneIDs_length,z=DBA_length)))


blasDBrRNA=blastFilterTranscriptType(blastDB,'rRNA',22)
rRNAIDs=queryIDs(blasDBrRNA)
length(rRNAIDs)
length(which(rRNAIDs %in% orfIdA))
length(which(rRNAIDs %in% orfIdB))
length(which(rRNAIDs %in% orfIdC))

DBA_rRNAIDs_length = unlist(lapply(lapply(DBA,function(x,y){which(x %in% y)}, y=rRNAIDs),length))
DBA_length=unlist(lapply(DBA,length))
length(unlist(lapply(which(DBA_rRNAIDs_length>0),function(x,y,z){if(y[x]==z[x]){x}},y=DBA_rRNAIDs_length,z=DBA_length)))



blasDBsnRNA=blastFilterTranscriptType(blastDB,'snRNA',22)
snRNAIDs=queryIDs(blasDBsnRNA)
length(snRNAIDs)
length(which(snRNAIDs %in% orfIdA))
length(which(snRNAIDs %in% orfIdB))
length(which(snRNAIDs %in% orfIdC))

blasDBsense_intronic=blastFilterTranscriptType(blastDB,'sense_intronic',22)
sense_intronicIDs=queryIDs(blasDBsense_intronic)
length(sense_intronicIDs)
length(which(sense_intronicIDs %in% orfIdA))
length(which(sense_intronicIDs %in% orfIdB))
length(which(sense_intronicIDs %in% orfIdC))

blasDB3primeNCRNA=blastFilterTranscriptType(blastDB,'3prime_overlapping_ncrna',22)
threeprimeNCRNAIDs=queryIDs(blasDB3primeNCRNA)
length(threeprimeNCRNAIDs)
length(which(threeprimeNCRNAIDs %in% orfIdA))
length(which(threeprimeNCRNAIDs %in% orfIdB))
length(which(threeprimeNCRNAIDs %in% orfIdC))


blasDBsense_overlapping=blastFilterTranscriptType(blastDB,'sense_overlapping',22)
sense_overlappingIDs=queryIDs(blasDBsense_overlapping)
length(sense_overlappingIDs)
length(which(sense_overlappingIDs %in% orfIdA))
length(which(sense_overlappingIDs %in% orfIdB))
length(which(sense_overlappingIDs %in% orfIdC))

blasDBnon_stop_decay=blastFilterTranscriptType(blastDB,'non_stop_decay',22)
non_stop_decayIDs=queryIDs(blasDBnon_stop_decay)
length(non_stop_decayIDs)
length(which(non_stop_decayIDs %in% orfIdA))
length(which(non_stop_decayIDs %in% orfIdB))
length(which(non_stop_decayIDs %in% orfIdC))

blasDBtranscribed_processed_pseudogene=blastFilterTranscriptType(blastDB,'transcribed_processed_pseudogene',22)
transcribed_processed_pseudogeneIDs=queryIDs(blasDBtranscribed_processed_pseudogene)
length(transcribed_processed_pseudogeneIDs)
length(which(transcribed_processed_pseudogeneIDs %in% orfIdA))
length(which(transcribed_processed_pseudogeneIDs %in% orfIdB))
length(which(transcribed_processed_pseudogeneIDs %in% orfIdC))

blasDBsnoRNA=blastFilterTranscriptType(blastDB,'snoRNA',22)
snoRNAIDs=queryIDs(blasDBsnoRNA)
length(snoRNAIDs)
length(which(snoRNAIDs %in% orfIdA))
length(which(snoRNAIDs %in% orfIdB))
length(which(snoRNAIDs %in% orfIdC))

blasDBunitary_pseudogene=blastFilterTranscriptType(blastDB,'unitary_pseudogene',22)
unitary_pseudogeneIDs=queryIDs(blasDBunitary_pseudogene)
length(unitary_pseudogeneIDs)
length(which(unitary_pseudogeneIDs %in% orfIdA))
length(which(unitary_pseudogeneIDs %in% orfIdB))
length(which(unitary_pseudogeneIDs %in% orfIdC))

blasDBmiRNA=blastFilterTranscriptType(blastDB,'miRNA',22)
miRNAIDs=queryIDs(blasDBmiRNA)
length(miRNAIDs)
length(which(miRNAIDs %in% orfIdA))
length(which(miRNAIDs %in% orfIdB))
length(which(miRNAIDs %in% orfIdC))

blasDBmiRNA=blastFilterTranscriptType(blastDB,'miRNA',22)
miRNAIDs=queryIDs(blasDBmiRNA)
length(miRNAIDs)
length(which(miRNAIDs %in% orfIdA))
length(which(miRNAIDs %in% orfIdB))
length(which(miRNAIDs %in% orfIdC))


blasDBpseudogene=blastFilterTranscriptType(blastDB,'pseudogene',22)
pseudogeneIDs=queryIDs(blasDBpseudogene)
length(pseudogeneIDs)
length(which(pseudogeneIDs %in% orfIdA))
length(which(pseudogeneIDs %in% orfIdB))
length(which(pseudogeneIDs %in% orfIdC))

blasDBTEC=blastFilterTranscriptType(blastDB,'TEC',22)
TECRNAIDs=queryIDs(blasDBTEC)
length(TECRNAIDs)
length(which(TECRNAIDs %in% orfIdA))
length(which(TECRNAIDs %in% orfIdB))
length(which(TECRNAIDs %in% orfIdC))




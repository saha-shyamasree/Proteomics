folder="C:/Data/HUMAN/mgf/"
folder="G:/SearchEngine_output/MSGF+/"
folder="C:/Users/shyama/Dropbox/results/"
folder="C:/Users/shyama/Dropbox/results/Bat_human_hendra/Human/"
filename_ut1="DM_from_raw_trinity_uniprot+fdr+th+grouping.csv"
msgf_ut1 = read.csv(file=paste(folder,filename_ut1,sep=""),header=TRUE)
filename_ut2="DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping.csv"
msgf_ut2 = read.csv(file=paste(folder,filename_ut2,sep=""),header=TRUE)
filename_t="DM_from_raw_trinity+fdr+th+grouping.csv"
filename_t="trinity_bat+fdr+th+grouping.csv"
msgf_t = read.csv(file=paste(folder,filename_t,sep=""),header=TRUE)
filename_u="DM_from_raw_uniprot+fdr+th+grouping.csv"
filename_u="uniprot_bat+fdr+th+grouping.csv"
msgf_u = read.csv(file=paste(folder,filename_u,sep=""),header=TRUE)

filename_ut1="DM_from_raw_trinity_uniprot+fdr+th+grouping+prt.csv"
prt_ut1 = read.csv(file=paste(folder,filename_ut1,sep=""),header=TRUE)
filename_ut2="DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping+prt.csv"
prt_ut2 = read.csv(file=paste(folder,filename_ut2,sep=""),header=TRUE)
#filename_t="DM_from_raw_trinity+fdr+th+grouping+prt.csv"
filename_t="trinity_bat+fdr+th+grouping+prt.csv"
prt_t = read.csv(file=paste(folder,filename_t,sep=""),header=TRUE)
prt_t_rev=prt_t[-(grep("_REVERSED",prt_t[,'protein.accession'])),]
prt_t_rev=as.matrix(prt_t_rev)

filename_u="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
filename_u="uniprot_bat+fdr+th+grouping+prt.csv"
prt_u = read.csv(file=paste(folder,filename_u,sep=""),header=TRUE)
prt_u_rev=prt_u[-(grep("_REVERSED",prt_u[,'protein.accession'])),]
prt_u_rev=as.matrix(prt_u_rev)
##########################  Protein Group / Proteins ################################
### unique peptide thresholding ###
threshold=1
prt_t_th = prt_t[which(prt_t$distinct.peptide.sequences>threshold),]
prt_t_th_rev=prt_t_th[-(grep("_REVERSED",prt_t_th[,'protein.accession'])),]
prt_t_th_rev=as.matrix(prt_t_th_rev)
prt_u_th = prt_u[which(prt_u$distinct.peptide.sequences>threshold),]
prt_u_th_rev=prt_u_th[-(grep("_REVERSED",prt_u_th[,'protein.accession'])),]
prt_u_th_rev=as.matrix(prt_u_th_rev)


prt_ut1_th = prt_ut1[which(prt_ut1$distinct.peptide.sequences>threshold),]
prt_ut2_th = prt_ut2[which(prt_ut2$distinct.peptide.sequences>threshold),]

############### protein ambiguity group ##########
length(unique(prt_t_th[,'PAG.ID']))
length(unique(prt_u_th[,'PAG.ID']))
length(unique(prt_ut1_th$PAG.ID))
length(unique(prt_ut2_th$PAG.ID))

######### Anchor protein count ##############
anc_prt_t = prt_t_th[grep('anchor protein',prt_t_th$group.membership),]
anc_prt_t=as.matrix(anc_prt_t)
anc_prt_t_rev=anc_prt_t[-(grep("_REVERSED",anc_prt_t[,'protein.accession'])),]
anc_prt_t_all = prt_t[grep('anchor protein',prt_u$group.membership),]
anc_prt_u = prt_u_th[grep('anchor protein',prt_u_th$group.membership),]
anc_prt_u=as.matrix(anc_prt_u)
anc_prt_u_rev=anc_prt_u[-(grep("_REVERSED",anc_prt_u[,'protein.accession'])),]

anc_prt_u_all = prt_u[grep('anchor protein',prt_u$group.membership),]
anc_prt_u_all_rev=anc_prt_u_all[-(grep("_REVERSED",anc_prt_u_all[,'protein.accession'])),]

mean(u_length[which(u_length[,1] %in% anc_prt_u_rev[,'protein.accession']),2])

anc_prt_ut1 = prt_ut1_th[grep('anchor protein',prt_ut1_th$group.membership),]
anc_prt_ut1=as.matrix(anc_prt_ut1)
anc_prt_ut2 = prt_ut2_th[grep('anchor protein',prt_ut2_th$group.membership),]
anc_prt_ut2=as.matrix(anc_prt_ut2)

length(which(anc_prt_t$protein.accession %in% duplicate_map$query_name))
length(which(anc_prt_u$description %in% duplicate_map$hit_def))
length(which(anc_prt_ut1$protein.accession %in% duplicate_map$query_name))
length(which(anc_prt_ut2$protein.accession %in% duplicate_map$query_name))

###### Protein Overlap calculation ######
d_folder="C:/Users/shyama/Dropbox/results/"
filename="duplicate_id_map_clean.csv"
duplicate_map = read.csv(file=paste(d_folder,filename,sep=""),header=TRUE)
### Uniprot ###
######Overlap with trinity####
uni_id = which(anc_prt_t[,'protein.accession'] %in% duplicate_map$query_name)
length(uni_d)
for(i in uni_id)
{
    anc_prt_t[i,'protein.accession'] = as.character(sub("([^\\s]+) (.*)?","\\1",as.character(duplicate_map[which(as.character(anc_prt_t[i, 'protein.accession']) == duplicate_map$query_name),'hit_def'])))
}

length(which(anc_prt_t[,'protein.accession'] %in% anc_prt_u[,'protein.accession']))  ##u and T
anc_prt_t_exc=anc_prt_t[which(!anc_prt_t[,'protein.accession'] %in% anc_prt_u[,'protein.accession']),]
anc_prt_u_exc=anc_prt_u[which(!anc_prt_u[,'protein.accession'] %in% anc_prt_t[,'protein.accession']),]

## Overlap, all protein, no threshold on peptide number ##

filename="blast_trinity_id_clean.csv"
filename="blastdb/blastCSV/p_alecto_blast.xml2.csv"
blast=read.csv(file=paste(folder,filename,sep=""),header=TRUE)
blast_map=subset(blast,match=='yes')
blast_res_t_exc=blast[which(blast_map$query_name %in% anc_prt_t_exc[,'protein.accession']),]
subset(blast_res_t_exc, $length_ratio==1 & $good_match>0.9 & $long_match>0.9)
blast_map_eval=subset(blast_map,e.value>=-0.9 & e.value<=0.1)
blast_map_eval=as.matrix(blast_map_eval)
blast_map_eval[,'hit_def']=sub("([^\\s]+) (.*)?","\\1",blast_map_eval[,'hit_def'])
blast_map_eval[,'hit_def']=sub("^sp","sw",blast_map_eval[,'hit_def'])


anc_prt_t_all_rev=anc_prt_t_all[-(grep("_REVERSED",anc_prt_t_all[,'protein.accession'])),]
anc_prt_t_all_rev = as.matrix(anc_prt_t_all_rev)
uni_id = which(anc_prt_t_all_rev[,'protein.accession'] %in% blast_map_eval[,'query_name'])
length(uni_id)
for(i in uni_id)
{
    anc_prt_t_all_rev[i,'protein.accession'] = blast_map_eval[which(blast_map_eval[,'query_name'] %in% anc_prt_t_all_rev[i,'protein.accession']),'hit_def']
}

length(which(anc_prt_t_all_rev[,'protein.accession'] %in% anc_prt_u_all_rev[,'protein.accession']))
anc_prt_t_all_rev_excl=anc_prt_t_all_rev[which(!anc_prt_t_all_rev[,'protein.accession'] %in% anc_prt_u_all_rev[,'protein.accession']),]

##proteins that are similar to that of uniprot proteins, but was only identified in Trinity
length(grep("^sw|^tr",anc_prt_t_all_rev_excl[,'protein.accession']))

head(anc_prt_t_all_rev_excl[grep("^sw|^tr",anc_prt_t_all_rev_excl[,'protein.accession']),])

## Overlap, all protein, threshold on peptide number>1 ##

anc_prt_t_rev=anc_prt_t[-(grep("_REVERSED",anc_prt_t[,'protein.accession'])),]
anc_prt_t_rev = as.matrix(anc_prt_t_rev)
anc_prt_t_all_rev[,'protein.accession'] = sub("([^,]+) (.*)?","\\1",anc_prt_t_all_rev[,'protein.accession'])
anc_prt_t_rev_b=anc_prt_t[-(grep("_REVERSED",anc_prt_t[,'protein.accession'])),]
anc_prt_t_rev_b = as.matrix(anc_prt_t_rev_b)

uni_id = which(anc_prt_t_rev[,'protein.accession'] %in% blast_map_eval[,'query_name'])

length(uni_id)
for(i in uni_id)
{
    anc_prt_t_rev[i,'protein.accession'] = blast_map_eval[which(blast_map_eval[,'query_name'] %in% anc_prt_t_rev[i,'protein.accession']),'hit_def']
}
length(which(sub("^(?:sw|tr)\\|(.*)?\\|(.*)?","\\1",unique(anc_prt_t_rev[,'protein.accession'])) %in% anc_prt_u_rev[,'protein.accession']))
t_u_ids=sub("^(?:sw|tr)\\|(.*)?\\|(.*)?","\\1",anc_prt_t_rev[,'protein.accession'])
anc_prt_t_rev_excl=anc_prt_t_rev[which(!sub("^(?:sw|tr)\\|(.*)?\\|(.*)?","\\1",anc_prt_t_rev[,'protein.accession']) %in% anc_prt_u_rev[,'protein.accession']),]

####bat
uni_id = which(anc_prt_t_rev[,'protein.accession'] %in% blast_map_eval[,'hit_def'])

length(uni_id)
for(i in uni_id)
{
    anc_prt_t_rev[i,'protein.accession'] = sub("([^ ]*)? (.*)","\\1",blast_map_eval[which(blast_map_eval[,'hit_def'] %in% anc_prt_t_rev[i,'protein.accession']),'query_name'])
}

blast_map_eval_qname=blast_map_eval
blast_map_eval_qname[,'query_name']=sub("([^ ]*)? (.*)","\\1",blast_map_eval[,'query_name'])
uni_id = which(anc_prt_u_rev[,'protein.accession'] %in% blast_map_eval_qname[,'query_name'])
prt_t_th_rev[,'protein.accession']=sub("([^,]*)?,(.*)","\\1",prt_t_th_rev[,'protein.accession'])
prt_t_rev[,'protein.accession']=sub("([^,]*)?,(.*)","\\1",prt_t_rev[,'protein.accession'])
length(uni_id)
for(i in uni_id)
{
    anc_prt_u_rev[i,'protein.accession'] = blast_map_eval_qname[which(blast_map_eval_qname[,'query_name'] %in% anc_prt_u_rev[i,'protein.accession']),'hit_def']
}
length(which(anc_prt_t_rev[,'protein.accession'] %in% anc_prt_u_rev[,'protein.accession']))

length(which(anc_prt_t_rev[,'protein.accession'] %in% anc_prt_u_rev[,'protein.accession']))
[1] 2548
> length(which(anc_prt_u_rev[,'protein.accession'] %in% anc_prt_t_rev[,'protein.accession']))
[1] 2683
> length(which(unique(anc_prt_u_rev[,'protein.accession']) %in% anc_prt_t_rev[,'protein.accession']))
[1] 2548
> length(which(unique(anc_prt_u_rev[,'protein.accession']) %in% prt_t_rev[,'protein.accession']))
[1] 3089
> length(which(unique(anc_prt_u_rev[,'protein.accession']) %in% prt_t_th_rev[,'protein.accession']))
[1] 3747


##proteins that are similar to that of uniprot proteins, but was only identified in Trinity
length(grep("^sw|^tr",unique(anc_prt_t_rev_excl[,'protein.accession'])))
anc_prt_t_rev_excl_uni=unique(anc_prt_t_rev_excl[grep("^sw|^tr",anc_prt_t_rev_excl[,'protein.accession']),'protein.accession'])
length(which(anc_prt_t_rev_excl_uni %in% prt_u_th_rev[,'protein.accession']))
length(which(anc_prt_t_rev_excl_uni %in% prt_u_all_rev[,'protein.accession']))
head(anc_prt_t_all_rev_excl[grep("^sw|^tr",anc_prt_t_all_rev_excl[,'protein.accession']),])
length(which(!unique(anc_prt_t_all_rev_excl[(grep("^sw|^tr",anc_prt_t_all_rev_excl[,'protein.accession'])),'protein.accession']) %in% prt_u_rev[,'protein.accession']))

##exclusively in uniprot

anc_prt_u_rev_excl=anc_prt_u_rev[which(!anc_prt_u_rev[,'protein.accession'] %in% unique(t_u_ids)),]
dim(anc_prt_u_rev_excl)

blast_map_eval[which(sub("^(?:sw|tr)\\|(.*)?\\|(.*)?","\\1",blast_map_eval[,'hit_def']) %in% 
uni_id = which(prt_t_th_rev[,'protein.accession'] %in% blast_map_eval[,'query_name'])
length(uni_id)
for(i in uni_id)
{
    prt_t_th_rev[i,'protein.accession'] = blast_map_eval[which(blast_map_eval[,'query_name'] %in% prt_t_th_rev[i,'protein.accession']),'hit_def']
}
prt_t_th_rev_unq_id=sub("^(?:sw|tr)\\|(.*)?\\|(.*)?","\\1",unique(prt_t_th_rev[,'protein.accession']))

prt_t_th_rev_unq_id_no_iso=sub("([^-]*)?-(.*)","\\1",prt_t_th_rev_unq_id)
anc_u_th_rev_no_iso=sub("([^-]*)?-(.*)?","\\1",unique(anc_prt_u_rev[,'protein.accession']))
length(which(anc_prt_u_rev_excl[,'protein.accession'] %in% prt_t_th_rev_unq_id))

for(i in 1:length(prt_t_th_rev_unq_id))
{
    grep(sub("([^-]*)","\\1",prt_t_th_rev_unq_id[i]),)
}





####### considering all protein with atleast 2 peptide

prt_u_th_rev[,'protein.accession']

## No filter ##
anc_t = prt_t[grep('anchor protein',prt_t$group.membership),]
anc_t=as.matrix(anc_t)
uni_id = which(anc_t[,'protein.accession'] %in% duplicate_map$query_name)
length(uni_id)
for(i in uni_id)
{
    anc_t[i,'protein.accession'] = as.character(sub("([^\\s]+) (.*)?","\\1",as.character(duplicate_map[which(as.character(anc_t[i, 'protein.accession']) == duplicate_map$query_name),'hit_def'])))
}

blast_res_t_exc_eval_hit=blast_res_t_exc_eval
blast_res_t_exc_eval_hit[,'hit_def']=sub("([^\\s]+) (.*)?","\\1",blast_res_t_exc_eval[,'hit_def'])


anc_u = prt_u[grep('anchor protein',prt_u$group.membership),]

length(which(anc_t[,'protein.accession'] %in% anc_u[,'protein.accession']))  ##u and T

####### Overlap with UT1 #######
uni_id = which(anc_prt_ut1[,'protein.accession'] %in% duplicate_map$query_name)
length(uni_d)
for(i in uni_id)
{
    anc_prt_ut1[i,'protein.accession'] = as.character(sub("([^\\s]+) (.*)?","\\1",as.character(duplicate_map[which(as.character(anc_prt_ut1[i, 'protein.accession']) == duplicate_map$query_name),'hit_def'])))
}
length(which(anc_prt_ut1[,'protein.accession'] %in% anc_prt_u[,'protein.accession']))  ##u and UT1

####### Overlap with UT2 #######
uni_id = which(anc_prt_ut2[,'protein.accession'] %in% duplicate_map$query_name)
length(uni_d)
for(i in uni_id)
{
    anc_prt_ut2[i,'protein.accession'] = as.character(sub("([^\\s]+) (.*)?","\\1",as.character(duplicate_map[which(as.character(anc_prt_ut2[i, 'protein.accession']) == duplicate_map$query_name),'hit_def'])))
}
length(which(anc_prt_ut2[,'protein.accession'] %in% anc_prt_u[,'protein.accession']))  ##u and UT2

### Trinity ###
length(which(anc_prt_ut1[,'protein.accession'] %in% anc_prt_t[,'protein.accession']))  ##t and UT1
length(which(anc_prt_ut2[,'protein.accession'] %in% anc_prt_t[,'protein.accession']))  ##t and UT2


##################### Peptide #############################
trinity_pep_seq = as.matrix(unique(msgf_t$Sequence))
uniprot_pep_seq = as.matrix(unique(msgf_u$Sequence))
ut1_pep_seq = as.matrix(unique(msgf_ut1$Sequence))
ut2_pep_seq = as.matrix(unique(msgf_ut2$Sequence))

length(trinity_pep_seq)
length(uniprot_pep_seq)
length(ut1_pep_seq)
length(ut2_pep_seq)

length(which(trinity_pep_seq %in% uniprot_pep_seq))
length(which(ut1_pep_seq %in% ut2_pep_seq))
length(which(ut2_pep_seq %in% ut1_pep_seq))
length(which(uniprot_pep_seq %in% trinity_pep_seq))
length(which(uniprot_pep_seq %in% ut1_pep_seq))
length(which(uniprot_pep_seq %in% ut2_pep_seq))
length(which(trinity_pep_seq %in% ut1_pep_seq))
length(which(trinity_pep_seq %in% ut2_pep_seq))


############################################## MaxQuant ########################################################
filename="peptides_trinity.txt"
mq_t=read.table(file=paste(folder,filename,sep=""),header=TRUE,sep="\t")
mq_t=as.matrix(mq_t)
blast_map=subset(blast,match=='yes')
inds=c()
blast_single_id = blast_map
blast_single_id$query_name = sub("^([^;]+)(.*)?","\\1",blast$query_name)
blast_single_id_th=subset(blast_single_id,e.value<=0)
inds=which(mq_t[,'Leading.razor.protein'] %in% blast_single_id_th$query_name)

for(i in inds)
{
    mq_t[i,'Leading.razor.protein'] = sub("sp\\|(.*)?\\|(.*)?|tr\\|(.*)?\\|(.*)?","\\1",blast_single_id_th[which(blast_single_id_th[,'query_name'] == mq_t[i,'Leading.razor.protein']),'hit_def'])       
}
filename="peptides_uniprot.txt"
mq_u=read.table(file=paste(folder,filename,sep=""),header=TRUE,sep="\t")
mq_u=as.matrix(mq_u)

mq_t=mq_t[-(grep("^CON_|^REV_",mq_t[,'Leading.razor.protein'])),]
mq_u=mq_u[-(grep("^CON_|^REV_",mq_u[,'Leading.razor.protein'])),]

unq_mq_u = unique(mq_u[,'Leading.razor.protein'])
unq_mq_t = unique(mq_t[,'Leading.razor.protein'])

u_length = read.table(file=paste(folder,"uniprot_id_length.tsv",sep=""),header=FALSE,sep="\t")
t_length = read.table(file=paste(folder,"trinity_id_length.tsv",sep=""),header=FALSE,sep="\t")
t_length_2=t_length
t_length_2[,1]=sub("^([^;]+)(.*)?","\\1",t_length_2[,1])
inds=which(t_length_2[,1] %in% blast_single_id_th$query_name)

for(i in inds)
{
    t_length_2[i,1] = sub("sp|tr\\|(.*)?\\|(.*)?","\\1",blast_single_id_th[which(blast_single_id_th[,'query_name'] == t_length_2[i,1]),'hit_def'])       
}
t_u_length=t_length
t_u_length=as.matrix(t_u_length)
uni_id = which(t_u_length[,1] %in% anc_prt_t_all_rev_backup[,'protein.accession'])
length(uni_id)
for(i in uni_id)
{
    ind=which(blast_map_eval[,'query_name'] %in% t_u_length[i,1])
    if(length(ind)>0)
    {
        t_u_length[i,1] = blast_map_eval[ind,'hit_def']   
    }
}

mean(t_length[which(t_length[,1] %in% anc_prt_t_all_rev_backup[,'protein.accession']),2])

u_length_ident=u_length[which(u_length[,1] %in% mq_u[,'Leading.razor.protein']),]
u_length_ident_msgf=u_length[which(u_length[,1] %in% [,'Leading.razor.protein']),]

which(u_length[,1] %in% sub("sw|tr\\|(.*)?\\|(.*)?","\\1",anc_prt_t_all_rev[,'protein.accession']))
temp=anc_prt_t_all_rev[grep("^sw|^tr",anc_prt_t_all_rev[,'protein.accession']),]
temp_op=anc_prt_t_all_rev[grep("^comp",anc_prt_t_all_rev[,'protein.accession']),]
unq_anc_prt_t_all_rev=unique(sub("^(?:sw|tr)\\|(.*)?\\|(.*)?","\\1",temp[,'protein.accession']))


(sum(t_length[which(t_length[,1] %in% unique(temp_op[,'protein.accession'])),2])+sum(u_length[which(u_length[,1] %in% unq_anc_prt_t_all_rev),2]))/(length(t_length[which(t_length[,1] %in% unique(temp_op[,'protein.accession'])),2])+length(u_length[which(u_length[,1] %in% unq_anc_prt_t_all_rev),2]))
length(unq_anc_prt_t_all_rev)

##################  Protein Grp MaxQuant###################
mq_prt_u = read.table(file=paste(folder,"proteinGroups_uniprot.txt",sep=""),header=TRUE,sep="\t")
mq_prt_u_th=subset(mq_prt_u,(Peptides>1))

mq_prt_u_th=mq_prt_u_th[-(c(2963,2964,2965,2966,2967,2968,2969,2970,2971,2972,2973)),]  ## these are REV protein grp. 411 has mostly target protein, but contains 1/2 REV.

filenameUAdeno="DM_from_raw_uniprot+fdr+th+grouping+prt.csv"
filenameTAdeno="trinity_PITORF+fdr+th+grouping+prt.csv"
filenameHUHendra="uniprot_human+fdr+th+grouping+prt.csv"
filenameHTHendra="trinity+fdr+th+grouping+prt.csv"
folderHHendra="C:/Users/shyama/Dropbox/results/Bat_human_hendra/Human/"
folderAdeno="C:/Users/shyama/Dropbox/results/Human_adenovirus/"

fbHHendra="uniprot_db_blast_human_hendra2.csv"
fbHfolder="C:/Users/shyama/Dropbox/results/blastdb/blastCSV/"
fbAdeno="trinity_PITORF_human_adeno_blast2.csv"
fbAdenofolder="C:/Users/shyama/Dropbox/results/blastdb/blastCSV/"

rev=1
peptide=1
pepThreshold=1
upper=0.1
lower=-0.9

UAdeno=proteinGroupFiltered(proteinGroup(filenameUAdeno,folderAdeno),rev,peptide,pepThreshold)
TAdeno=proteinGroupFiltered(proteinGroup(filenameTAdeno,folderAdeno),rev,peptide,pepThreshold)
TAdeno=replaceComma(TAdeno,3)
bAdeno=blastFilterEval(blastFilterMatch(blast(fbAdeno,fbAdenofolder),1),upper,lower)
bAdeno=firstWord(bAdeno,4)


HUHendra=proteinGroupFiltered(proteinGroup(filenameHUHendra,folderHHendra),rev,peptide,pepThreshold)
HTHendra=proteinGroupFiltered(proteinGroup(filenameHTHendra,folderHHendra),rev,peptide,pepThreshold)
HTHendra=replaceComma(HTHendra,3)
bHHendra=blastFilterEval(blastFilterMatch(blast(fbHHendra,fbHfolder),1),upper,lower)
bHHendra=firstWord(bHHendra,4)

length(overlap(proteinGrpAnchor(UAdeno)[,'protein.accession'],HUHendra))
dim(proteinGrpAnchor(UAdeno))
length(overlap(proteinGrpAnchor(HUHendra)[,'protein.accession'],UAdeno))
dim(proteinGrpAnchor(HUHendra))

length(grep("OS=Hendra",proteinGrpAnchor(HUHendra)[,'description']))
UAdeno[,'protein.accession']=sub("^sw","sp",UAdeno[,'protein.accession'])
HUHendra[,'protein.accession']=sub("^sw","sp",HUHendra[,'protein.accession'])
HTHendra[,'protein.accession']=replaceIds(HTHendra,bHHendra)
TAdeno[,'protein.accession']=replaceIds(TAdeno,bAdeno)
length(unique(proteinGrpAnchor(HTHendra)[,'protein.accession']))
length(unique(proteinGrpAnchor(TAdeno)[,'protein.accession']))
length(proteinGrpAnchor(TAdeno)[,'protein.accession'])
length(overlap(proteinGrpAnchor(UAdeno)[,'protein.accession'],HUHendra))
length(overlap(proteinGrpAnchor(TAdeno)[,'protein.accession'],HTHendra))
length(overlap(proteinGrpAnchor(HTHendra)[,'protein.accession'],TAdeno))
length(overlap(proteinGrpAnchor(HTHendra)[,'protein.accession'],HUHendra))
length(overlap(proteinGrpAnchor(HUHendra)[,'protein.accession'],HTHendra))
length(grep("OS=Hendra",proteinGrpAnchor(HUHendra)[,'description']))
length(grep("OS=Hendra",proteinGrpAnchor(HTHendra)[,'description']))
length(grep("OS=Hendra",proteinGrpAnchor(HUHendra)[overlap(proteinGrpAnchor(HUHendra)[,'protein.accession'],HTHendra),'description']))
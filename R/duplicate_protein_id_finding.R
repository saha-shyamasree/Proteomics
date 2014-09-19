duplicate_trinity = read.csv("E:/Data/HUMAN/database/Trinity/duplicate_trinity_ids.csv",header=TRUE)
m=gregexpr("comp\\S+", duplicate_trinity$query_name, perl=TRUE)
duplicate_trinity_ids=unlist(regmatches(duplicate_trinity$query_name,m))


#########################################################################
pep_t=read.delim(file="E:/Data/HUMAN/MS/Results/Peptide_identification/MaxQuant/peptides_trinity.txt",header=TRUE)
pep_u=read.delim(file="E:/Data/HUMAN/MS/Results/Peptide_identification/MaxQuant/peptides_uniprot.txt",header=TRUE)
pep_ut1=read.delim(file="E:/Data/HUMAN/MS/Results/Peptide_identification/MaxQuant/peptides_uniprot_trinity.txt",header=TRUE)
pep_ut2=read.delim(file="E:/Data/HUMAN/MS/Results/Peptide_identification/MaxQuant/peptides_uniprot_trinity_no_duplicate.txt",header=TRUE)

pep_ut12_1=pep_ut1[which(pep_ut1$Sequence %in% pep_ut2$Sequence),] #_1 is added because, columns other than 'Sequence' might have different value for ut2.
pep_ut12_2=pep_ut2[which(pep_ut2$Sequence %in% pep_ut1$Sequence),] #_2 is added because, columns other than 'Sequence' might have different value for ut1.
pep_ut12_2$Leading.razor.protein=sub(".*\\|(.*)\\|.*$","\\1",pep_ut12_2$Leading.razor.protein,perl=TRUE)
pep_ut12_1[2,c("Sequence","Mass","Proteins","Leading.razor.protein")] 

##########################################################################
duplicate_id_map=read.csv(file="E:/Data/HUMAN/database/Trinity/comparison/duplicate_id_map.csv",header=T)
duplicate_id_map$query_name=sub("([^\\s]*)(.*)$","\\1",duplicate_id_map$query_name,perl=TRUE)
duplicate_id_map$hit_def=sub(".*\\|(.*)\\|.*$","\\1",duplicate_id_map$hit_def,perl=TRUE)
ut12_1_seq_prot = pep_ut12_1[,c("Sequence","Leading.razor.protein")]
ut12_2_seq_prot = pep_ut12_2[,c("Sequence","Leading.razor.protein")]
res=list()
for(i in 1:dim(ut12_1_seq_prot)[1])
{
    res[[i]]=which(as.character(ut12_1_seq_prot[i,1])==ut12_2_seq_prot[,1] & as.character(ut12_1_seq_prot[i,2])==ut12_2_seq_prot[,2])
}
ut12_1_seq_prot = as.matrix(ut12_1_seq_prot)
ut12_2_seq_prot = as.matrix(ut12_2_seq_prot)
count=0
count1=0
count2=0
count3=0 #duplicate protein
count4=0
count5=0 #one uni and one trinity
indices_both_trinity=c()
indices_both_uniprot=c()
for(i in 1:length(res))
{
    if(length(res[[i]])==1)
    {
        count=count+1;
    }
    else
    {
        if(length(res[[i]])==0)
        {
            ##check whether id is from different DB, if starts with 'comp', its trinity id, then check whether it is one of those duplicate proteins, if so, what is uniprot mapping? does that match? if the id is from same db(either trinity or uniprot), then investigate why different protein was assigned?
            if(ut12_1_seq_prot[i,1]==ut12_2_seq_prot[i,1])
            {
                if(length(grep("^(?!comp).+",ut12_1_seq_prot[i,2],perl=TRUE))>0 & length(grep("^(?!comp).+",ut12_2_seq_prot[i,2],perl=TRUE))>0)
                {
                    #both uniprot, but different
                    print(paste("Both uniprot",ut12_1_seq_prot[i,2],ut12_2_seq_prot[i,2],sep=","))
                    count1=count1+1
                }
                else
                {
                    if(length(grep("^comp.+",ut12_1_seq_prot[i,2],perl=TRUE))>0 & length(grep("^comp.+",ut12_2_seq_prot[i,2],perl=TRUE))>0)
                    {
                        #both trinity, but different
                        print(paste("Both trinity",ut12_1_seq_prot[i,2],ut12_2_seq_prot[i,2],sep=","))
                        count2=count2+1
                    }
                    else
                    {
                        print(paste("one uni and other trinity",ut12_1_seq_prot[i,2],ut12_2_seq_prot[i,2],sep=","))
                        count5=count5+1
                        if(length(grep("^comp",ut12_1_seq_prot[i,2],perl=TRUE))>0)
                        {
                            index=which(duplicate_id_map$query_name==ut12_1_seq_prot[i,2])
                            if(length(duplicate_id_map[index,"hit_def"] == ut12_2_seq_prot[i,2])>0)
                            {
                                count3=count3+1
                            }
                            else
                            {
                                count4=count4+1
                            }
                        }
                        else
                        {
                            if(length(grep("^comp",ut12_2_seq_prot[i,2],perl=TRUE)))
                            {
                                index=which(duplicate_id_map$query_name==ut12_2_seq_prot[i,2])
                                if(length(duplicate_id_map[index,"hit_def"] == ut12_1_seq_prot[i,2])>0)
                                {
                                    count3=count3+1
                                }
                                else
                                {
                                    count4=count4+1
                                }
                            }
                            else{
                                print("should be error")
                            }
                        }
                    }
                }
            }
        }
        else{
            #should not happen for peptide matching
            print("ERROR")
        }
    }
}
#res=apply(ut12_1_seq_prot,1,function(x){which(as.character(x[1])==ut12_2_seq_prot[,1] & as.character(x[2])==ut12_2_seq_prot[,2])})


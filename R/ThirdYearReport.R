##Paper Analysis of FPKM values
## Spearman Rank Test to find feature / column correlation 
	cor_val = cor( a, method = "spearman" )
	require(pheatmap)
	pheatmap( cor_val, display_numbers = T, cluster_rows= F, cluster_cols = F, fontsize_number = 12  )
	number_format

## Pair wise correlation plot per feature  
	panel.cor <- function(x, y, digits = 2, cex.cor, ...)
	{
        usr = par("usr"); on.exit(par(usr))
        par(usr=c(0,1,0,1))
        corcoff=cor.test(x,y, method = "spearman")
      text( 0.5, 0.4, paste( "Cor:", round(corcoff$estimate, digits=2), sep = "" ), cex = 1.4)

    }

#temp=read.csv(file="C:/Users/Shyamasree/Dropbox/PITDB/HumanAdeno/Score/human_adeno2.csv",stringsAsFactors=FALSE)
#fpkm=read.csv(file="G:/Bristol/Human/adenovirus/RSEM/human_adeno.isoforms.results.identified.tsv", sep="\t")

#temp=read.csv(file="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Score/human_adeno2.csv")
#temp_fpkm=cbind(temp,FPKM=rep_len(0,length(feats$PepQ_Value)))
temp_fpkm=temp
temp_fpkm$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp_fpkm$TGE_Id), '|', fixed=TRUE)))[,1]
temp_fpkm$FPKM<-fpkm[match(temp_fpkm$transcript_id, fpkm$transcript_id),'FPKM']

feats=temp_fpkm[c('PepQ_Value','PepCount','PSMCount','UniquePepCount','PeptideCoverage','FPKM')]
feats$PepQ_Value=as.numeric(as.character(feats$PepQ_Value))
feats$PepCount=as.numeric(as.character(feats$PepCount))
feats$PSMCount=as.numeric(as.character(feats$PSMCount))
feats$UniquePepCount=as.numeric(as.character(feats$UniquePepCount))
feats$PeptideCoverage=as.numeric(as.character(feats$PeptideCoverage))
feats$FPKM=as.numeric(as.character(feats$FPKM))

pairs(feats, upper.panel = panel.cor,cex.label=1.5)

cor_val = cor( feats, method = "spearman" )
require(pheatmap)
pheatmap( cor_val, display_numbers = T, cluster_rows= F, cluster_cols = F, fontsize_number = 12  )
number_format

temp_iso=temp_fpkm[grepl("_",temp$TGE_Class),]
feats_iso=temp_iso[c('PepQ_Value','PepCount','PSMCount','UniquePepCount','PeptideCoverage','FPKM')]
feats_iso$PepQ_Value=as.numeric(as.character(feats_iso$PepQ_Value))
feats_iso$PepCount=as.numeric(as.character(feats_iso$PepCount))
feats_iso$PSMCount=as.numeric(as.character(feats_iso$PSMCount))
feats_iso$UniquePepCount=as.numeric(as.character(feats_iso$UniquePepCount))
feats_iso$PeptideCoverage=as.numeric(as.character(feats_iso$PeptideCoverage))
feats_iso$FPKM=as.numeric(as.character(feats_iso$FPKM))
pairs(feats_iso, upper.panel = panel.cor,cex.label=1.5)

library(sm)

temp_fpkm$FPKM=as.numeric(as.character(temp_fpkm$FPKM))
#temp$Score=as.numeric(as.character(temp$Score))
temp_fpkm$TGE_Class=as.character(temp_fpkm$TGE_Class)
tgeClass=temp_fpkm$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isoform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)

# plot densities
sm.density.compare(temp_fpkm$FPKM, tgeClass.f, xlab="FPKM")
title(main="Human Adeno")

library(ggplot2)
temp_fpkm_backup=temp_fpkm
temp_fpkm[temp_fpkm$FPKM>=100,'FPKM']=100
# First type of color
p <- ggplot(temp_fpkm, aes(factor(tgeClass), FPKM))
p + geom_violin(aes(fill = tgeClass))
dev.off()
 
# Second type
p <- ggplot(temp_fpkm, aes(factor(tgeClass), FPKM))
p + geom_violin(aes(fill = factor(tgeClass)))

################################################################################################################
##Third Year report
a = rnorm(1000)
dim(a) = c(200,5)

## Data set with 200 observations and 5 feature (columns) ##

## Spearman Rank Test to find feature / column correlation 
	cor_val = cor( a, method = "spearman" )
	require(pheatmap)
	pheatmap( cor_val, display_numbers = T, cluster_rows= F, cluster_cols = F, fontsize_number = 12  )
	number_format

## Pair wise correlation plot per feature  
	panel.cor <- function(x, y, digits = 2, cex.cor, ...)
	{
        usr = par("usr"); on.exit(par(usr))
        par(usr=c(0,1,0,1))
        corcoff=cor.test(x,y, method = "spearman")
      text( 0.5, 0.4, paste( "Cor:", round(corcoff$estimate, digits=2), sep = "" ), cex = 1.4)

    }
	pairs(a, upper.panel = panel.cor,cex.label=1.5)
temp[,3]=as.numeric(temp[,3])
temp[,4]=as.numeric(temp[,4])
temp[,5]=as.numeric(temp[,5])
temp[,6]=as.numeric(temp[,6])
temp[,7]=as.numeric(temp[,7])
apply(temp[,c(3,4,5,6,7)], 1, function(x) sum(x) )

scale<-function(data, bottom = 0, top = 20){
    min1 = min(data);
    max1 = max(data);
    scaled_data = ( data - min1 ) *  ( ( top - bottom ) / ( max1- min1 ) ) +  bottom;
    return(scaled_data);
} 


#temp=read.csv(file="C:/Users/shyama/Dropbox/PITDB/HumanAdeno/Score/human_adeno2.csv")
temp=read.csv(file="C:/Users/Shyamasree/Dropbox/PITDB/HumanAdeno/Score/human_adeno2.csv",stringsAsFactors=FALSE)
fpkm=read.csv(file="G:/Bristol/Human/adenovirus/RSEM/human_adeno.isoforms.results.identified.tsv", sep="\t")

#temp=read.csv(file="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Score/human_adeno2.csv")
feats=temp[c('PepQ_Value','PepCount','PSMCount','UniquePepCount','PeptideCoverage')]
feats$PepQ_Value=as.numeric(as.character(feats$PepQ_Value))
feats$PepCount=as.numeric(as.character(feats$PepCount))
feats$PSMCount=as.numeric(as.character(feats$PSMCount))
feats$UniquePepCount=as.numeric(as.character(feats$UniquePepCount))
feats$PeptideCoverage=as.numeric(as.character(feats$PeptideCoverage))
feats=cbind(feats,FPKM=rep_len(0,length(feats$PepQ_Value)))

pairs(feats, upper.panel = panel.cor,cex.label=1.5)

temp_iso=temp[grepl("_",temp$TGE_Class),]
feats_iso=temp_iso[c('PepQ_Value','PepCount','PSMCount','UniquePepCount','PeptideCoverage')]
feats_iso$PepQ_Value=as.numeric(as.character(feats_iso$PepQ_Value))
feats_iso$PepCount=as.numeric(as.character(feats_iso$PepCount))
feats_iso$PSMCount=as.numeric(as.character(feats_iso$PSMCount))
feats_iso$UniquePepCount=as.numeric(as.character(feats_iso$UniquePepCount))
feats_iso$PeptideCoverage=as.numeric(as.character(feats_iso$PeptideCoverage))
temp_fpkm=cbind(temp,FPKM=rep_len(0,length(feats$PepQ_Value)))

temp_fpkm$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp_fpkm$TGE_Id), '|', fixed=TRUE)))[,1]
temp_fpkm$FPKM<-fpkm[match(temp_fpkm$transcript_id, fpkm$transcript_id),'FPKM']
pairs(feats_iso, upper.panel = panel.cor,cex.label=1.5)


######################Frequency distribution of q-value #####################
library(sm)
#temp=read.csv(file="C:/Users/shyama/Dropbox/PITDB/HumanAdeno/Score/human_adeno2.csv")
temp=read.csv(file="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/Score/human_adeno2.csv")
temp$PepQ_Value=1-as.numeric(as.character(temp$PepQ_Value))
temp$Score=as.numeric(as.character(temp$Score))
temp$TGE_Class=as.character(temp$TGE_Class)
tgeClass=temp$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)

# plot densities
sm.density.compare(temp$PepQ_Value, tgeClass.f, xlab="Miles Per Gallon")
title(main="Human Adeno")

#Score
sm.density.compare(temp$Score, tgeClass.f, xlab="Score")
title(main="Human Adeno")


# add legend via mouse click
colfill<-c(2:(2+length(levels(tgeClass.f))))
legend(locator(1), levels(tgeClass.f), fill=colfill)


#######All isoform class###########
tgeClass.f <- factor(temp$TGE_Class)

# plot densities
sm.density.compare(temp$PepQ_Value, tgeClass.f, xlab="PepQ-Value")
colfill<-c(2:(2+length(levels(tgeClass.f))))
legend(locator(1), levels(tgeClass.f), fill=colfill)


#Score
sm.density.compare(temp$Score, tgeClass.f, xlab="Score")
title(main="Human Adeno")
colfill<-c(2:(2+length(levels(tgeClass.f))))
legend(locator(1), levels(tgeClass.f), fill=colfill)

#############################################################################################################################################
        ### Mouse ###

library(sm)
#dir="C:/Users/shyama/Dropbox/PITDB/MouseNelsonBay/Score/"
#dir="C:/Users/Shyamasree/Dropbox/PITDB/MouseNelsonBay/Score/"
#######################################################

##Paper
dirS="G:/Bristol/Mouse/PITDB/Score/"
dirR="G:/Bristol/Mouse/RSEM/"
files=list.files(path=dirS,pattern="*_score.csv")
for(i in 1:length(files))
{
    sam=gsub("_score.csv","",files[i])
    if(i==1)
    {
        temp=read.csv(file=paste(dirS,files[i],sep=""), stringsAsFactors = F )
        fpkm=read.csv(file=paste(dirR,'/',sam,".isoforms.results.identified.tsv",sep=""), sep="\t", stringsAsFactors = F )
        temp$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp$TGE_Id), '|', fixed=TRUE)))[,1]
        temp$FPKM<-fpkm[match(temp$transcript_id, fpkm$transcript_id),'FPKM']
        temp_fpkm=temp
    }
    else
    {
        temp=read.csv(file=paste(dirS,files[i],sep=""), stringsAsFactors = F)
        fpkm=read.csv(file=paste(dirR,'/',sam,".isoforms.results.identified.tsv",sep=""), sep="\t", stringsAsFactors = F )
        temp$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp$TGE_Id), '|', fixed=TRUE)))[,1]
        temp$FPKM<-fpkm[match(temp$transcript_id, fpkm$transcript_id),'FPKM']
        temp_fpkm=rbind(temp_fpkm,temp)
    }
}
feats=temp_fpkm[c('PepQ_Value','PepCount','PSMCount','UniquePepCount','PeptideCoverage','FPKM')]
feats$PepQ_Value=as.numeric(as.character(feats$PepQ_Value))
feats$PepCount=as.numeric(as.character(feats$PepCount))
feats$PSMCount=as.numeric(as.character(feats$PSMCount))
feats$UniquePepCount=as.numeric(as.character(feats$UniquePepCount))
feats$PeptideCoverage=as.numeric(as.character(feats$PeptideCoverage))
feats$FPKM=as.numeric(as.character(feats$FPKM))

pairs(feats, upper.panel = panel.cor,cex.label=1.5)

temp_fpkm$FPKM=as.numeric(as.character(temp_fpkm$FPKM))
#temp$Score=as.numeric(as.character(temp$Score))
temp_fpkm$TGE_Class=as.character(temp_fpkm$TGE_Class)
tgeClass=temp_fpkm$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isoform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)
temp_fpkm_backup=temp_fpkm
temp_fpkm[temp_fpkm$FPKM>=100,'FPKM']=100
# First type of color
p <- ggplot(temp_fpkm, aes(factor(tgeClass), FPKM))
p + geom_violin(aes(fill = tgeClass))

####################################
temp$PepQ_Value=1-as.numeric(as.character(temp$PepQ_Value))
temp$Score=as.numeric(as.character(temp$Score))
temp$TGE_Class=as.character(temp$TGE_Class)
tgeClass=temp$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)
#Score
sm.density.compare(temp$Score, tgeClass.f, xlab="Score")
title(main="Mouse")
colfill<-c(2:(2+length(levels(tgeClass.f))))
legend(locator(1), levels(tgeClass.f), fill=colfill)

#Score
par(mfcol=c(2,3))
for(i in 3:dim(temp)[2])
{
    #plot(density(as.numeric(as.character(temp[,i]))))
    sm.density.compare(as.numeric(as.character(temp[,i])), tgeClass.f, xlab=colnames(temp)[i])
    #title(main="Mouse")
    colfill<-c(2:(2+length(levels(tgeClass.f))))
    legend(locator(1), levels(tgeClass.f), fill=colfill)
   
}

#############################################################################
##### BaT


library(sm)
dir="C:/Users/shyama/Dropbox/PITDB/BatNelsonBay/Score/"
files=list.files(path=dir,pattern="*_score.csv")
for(i in 1:length(files))
{
    if(i==1)
    {
        temp=read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F )
    }
    else
    {
        temp=rbind(temp,read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F))
    }
}

temp$PepQ_Value=1-as.numeric(as.character(temp$PepQ_Value))
temp$Score=as.numeric(as.character(temp$Score))
temp$TGE_Class=as.character(temp$TGE_Class)
tgeClass=temp$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)
#Score
#sm.density.compare(temp$Score, tgeClass.f, xlab="Score")
#title(main="Bat")
#colfill<-c(2:(2+length(levels(tgeClass.f))))
#legend(locator(1), levels(tgeClass.f), fill=colfill)

#Score
par(mfcol=c(2,3))
for(i in 3:dim(temp)[2])
{
    #plot(density(as.numeric(as.character(temp[,i]))))
    sm.density.compare(as.numeric(as.character(temp[,i])), tgeClass.f, xlab=colnames(temp)[i])
    title(main="Bat")
    colfill<-c(2:(2+length(levels(tgeClass.f))))
    legend(locator(1), levels(tgeClass.f), fill=colfill)
   
}


#############################################################################
##Paper, with FPKM

dirS="G:/Bristol/Bat/NelsonBay/PITDB/Score/"
dirR="G:/Bristol/Bat/NelsonBay/RSEM/"
files=list.files(path=dirS,pattern="*_score.csv")
for(i in 1:length(files))
{
    sam=gsub("_score.csv","",files[i])
    if(i==1)
    {
        temp=read.csv(file=paste(dirS,files[i],sep=""), stringsAsFactors = F )
        fpkm=read.csv(file=paste(dirR,'/',sam,".isoforms.results.identified.tsv",sep=""), sep="\t", stringsAsFactors = F )
        temp$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp$TGE_Id), '|', fixed=TRUE)))[,1]
        temp$FPKM<-fpkm[match(temp$transcript_id, fpkm$transcript_id),'FPKM']
        temp_fpkm=temp
    }
    else
    {
        temp=read.csv(file=paste(dirS,files[i],sep=""), stringsAsFactors = F)
        fpkm=read.csv(file=paste(dirR,'/',sam,".isoforms.results.identified.tsv",sep=""), sep="\t", stringsAsFactors = F )
        temp$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp$TGE_Id), '|', fixed=TRUE)))[,1]
        temp$FPKM<-fpkm[match(temp$transcript_id, fpkm$transcript_id),'FPKM']
        temp_fpkm=rbind(temp_fpkm,temp)
    }
}
feats=temp_fpkm[c('PepQ_Value','PepCount','PSMCount','UniquePepCount','PeptideCoverage','FPKM')]
feats$PepQ_Value=as.numeric(as.character(feats$PepQ_Value))
feats$PepCount=as.numeric(as.character(feats$PepCount))
feats$PSMCount=as.numeric(as.character(feats$PSMCount))
feats$UniquePepCount=as.numeric(as.character(feats$UniquePepCount))
feats$PeptideCoverage=as.numeric(as.character(feats$PeptideCoverage))
feats$FPKM=as.numeric(as.character(feats$FPKM))

pairs(feats, upper.panel = panel.cor,cex.label=1.5)

temp_fpkm$FPKM=as.numeric(as.character(temp_fpkm$FPKM))
#temp$Score=as.numeric(as.character(temp$Score))
temp_fpkm$TGE_Class=as.character(temp_fpkm$TGE_Class)
tgeClass=temp_fpkm$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isoform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)
temp_fpkm_backup=temp_fpkm
temp_fpkm[temp_fpkm$FPKM>=100,'FPKM']=100
# First type of color
p <- ggplot(temp_fpkm, aes(factor(tgeClass), FPKM))
p + geom_violin(aes(fill = tgeClass))

#############################################################################
#Oliver's data

dir="C:/Users/shyama/Dropbox/PITDB/OliverData/Score/"
files=list.files(path=dir,pattern="*_score.csv")
for(i in 1:length(files))
{
    if(i==1)
    {
        temp=read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F )
    }
    else
    {
        temp=rbind(temp,read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F))
    }
}

temp$PepQ_Value=1-as.numeric(as.character(temp$PepQ_Value))
temp$Score=as.numeric(as.character(temp$Score))
temp$TGE_Class=as.character(temp$TGE_Class)
tgeClass=temp$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)


#Score
par(mfcol=c(2,3))
for(i in 3:dim(temp)[2])
{
    #plot(density(as.numeric(as.character(temp[,i]))))
    sm.density.compare(as.numeric(as.character(temp[,i])), tgeClass.f, xlab=colnames(temp)[i])
    title(main="Human Ovarian Cancer")
    colfill<-c(2:(2+length(levels(tgeClass.f))))
    legend(locator(1), levels(tgeClass.f), fill=colfill)
   
}
########################################################################################
##Paper
dirS="G:/Oliver/PITDB/Score/"
dirR="G:/Oliver/RSEM/"
files=list.files(path=dirS,pattern="*_score.csv")
for(i in 1:length(files))
{
    sam=gsub("_score.csv","",files[i])
    if(file.exists(paste(dirR,sam,".isoforms.results.identified.tsv")))
    {
        if(i==1)
        {
            temp=read.csv(file=paste(dirS,files[i],sep=""), stringsAsFactors = F )
            fpkm=read.csv(file=paste(dirR,sam,".isoforms.results.identified.tsv",sep=""), sep="\t", stringsAsFactors = F )
            temp$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp$TGE_Id), '|', fixed=TRUE)))[,1]
            temp$FPKM<-fpkm[match(temp$transcript_id, fpkm$transcript_id),'FPKM']
            temp_fpkm=temp
        }
        else
        {
            temp=read.csv(file=paste(dirS,files[i],sep=""), stringsAsFactors = F)
            fpkm=read.csv(file=paste(dirR,sam,".isoforms.results.identified.tsv",sep=""), sep="\t", stringsAsFactors = F )
            temp$transcript_id=data.frame(do.call('rbind', strsplit(as.character(temp$TGE_Id), '|', fixed=TRUE)))[,1]
            temp$FPKM<-fpkm[match(temp$transcript_id, fpkm$transcript_id),'FPKM']
            temp_fpkm=rbind(temp_fpkm,temp)
        }
    }
}
feats=temp_fpkm[c('PepQ_Value','PepCount','PSMCount','UniquePepCount','PeptideCoverage','FPKM')]
feats$PepQ_Value=as.numeric(as.character(feats$PepQ_Value))
feats$PepCount=as.numeric(as.character(feats$PepCount))
feats$PSMCount=as.numeric(as.character(feats$PSMCount))
feats$UniquePepCount=as.numeric(as.character(feats$UniquePepCount))
feats$PeptideCoverage=as.numeric(as.character(feats$PeptideCoverage))
feats$FPKM=as.numeric(as.character(feats$FPKM))

pairs(feats, upper.panel = panel.cor,cex.label=1.5)

temp_fpkm$FPKM=as.numeric(as.character(temp_fpkm$FPKM))
#temp$Score=as.numeric(as.character(temp$Score))
temp_fpkm$TGE_Class=as.character(temp_fpkm$TGE_Class)
tgeClass=temp_fpkm$TGE_Class
tgeClass[grep("_",tgeClass)] <- "isoform"
#temp$TGE_Class[grepl("_",temp$TGE_Class)] <- "isoform"
#temp[is.na(temp)]="isoform"
# create value labels
tgeClass.f <- factor(tgeClass)
temp_fpkm_backup=temp_fpkm
temp_fpkm[temp_fpkm$FPKM>=100,'FPKM']=100
# First type of color
p <- ggplot(temp_fpkm, aes(factor(tgeClass), FPKM))
p + geom_violin(aes(fill = tgeClass))
########################################################################################
###Mosquito

temp=read.csv(file="C:/Users/shyama/Dropbox/PITDB/Mosquito/Score/human_adeno2.csv")

########################################################################################


###################### Staked barplot #######################################
data=read.csv(file="C:/Users/shyama/Dropbox/Report/ThirdYear/dataAnalysis/canonical_trembl_isoform.csv", row.names = 1)
counts=t(data[c(1,2,3)])
row.names(counts)=c("Reviewed Canonical","Tr-EMBL","Reviewed Isoform")
barplot(counts, main="Known protein distribution",xlab="Datasets", col=c("grey28","grey40","grey84"), legend = rownames(counts))

#############################################################################
#Bat Counting total number of Uniprot proteins.

library(sm)
dir="C:/Users/shyama/Dropbox/PITDB/BatNelsonBay/Annotation/"
files=list.files(path=dir,pattern="*_details_annotation.csv")
for(i in 1:length(files))
{
    if(i==1)
    {
        temp=read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F )
    }
    else
    {
        temp=rbind(temp,read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F))
    }
}

length(unique(temp[temp$Source=='tr' & temp$Class=="known",'Protein.ID']))
length(unique(temp[temp$Source=='sp' & temp$Class=="known",'Protein.ID']))
length(unique(temp[temp$Class=="known",'Protein.ID']))
length(unique(temp$Protein.ID))

#Mouse tr = 922
#Mouse sp = 1462
#Bat tr = 764
#Human tr=109
#HUman sp=1747
#Mosquito
########################## All data ###########################################
temp=read.csv(file="C:/Users/shyama/Dropbox/PITDB/HumanAdeno/Score/human_adeno2.csv", stringsAsFactors = F)
dir="C:/Users/shyama/Dropbox/PITDB/OliverData/Score/"
files=list.files(path=dir,pattern="*_score.csv")
for(i in 1:length(files))
{
    
    temp=rbind(temp,read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F))
}
dir="C:/Users/shyama/Dropbox/PITDB/BatNelsonBay/Score/"
files=list.files(path=dir,pattern="*_score.csv")
for(i in 1:length(files))
{
    temp=rbind(temp,read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F))
}
dir="C:/Users/shyama/Dropbox/PITDB/MouseNelsonBay/Score/"
files=list.files(path=dir,pattern="*_score.csv")
for(i in 1:length(files))
{
    temp=rbind(temp,read.csv(file=paste(dir,files[i],sep=""), stringsAsFactors = F))
}
temp=rbind(temp,read.csv(file="C:/Users/shyama/Dropbox/PITDB/Mosquito/Score/aedes_score.csv", stringsAsFactors = F))

temp[,3]=as.numeric(temp[,3])
temp[,4]=as.numeric(temp[,4])
temp[,5]=as.numeric(temp[,5])
temp[,6]=as.numeric(temp[,6])
temp[,7]=as.numeric(temp[,7])
temp_back_up=temp

#temp[,3]=scale(temp[,3],0,1)
temp[,4]=scale(temp[,4],0,1)
temp[,5]=scale(temp[,5],0,1)
temp[,6]=scale(temp[,6],0,1)
temp[,7]=scale(temp[,7],0,1)
temp$PepQ_Value=1-temp$PepQ_Value
Score=apply(temp[,c(3,4,5,6,7)], 1, function(x) sum(x) )
temp_score=cbind(temp,Score)
temp_new=temp_score[,c(1:7,9)]
temp=temp_new
#temp[,8]=colSums(temp[,3:7])
temp$PepQ_Value=1-as.numeric(temp$PepQ_Value)
tgeClass=temp$TGE_Class
temp$TGE_Class[grep("_",temp$TGE_Class)] <- "isform"

tgeClass.f <- factor(temp$TGE_Class)
par(mfcol=c(2,3))
for(i in 3:dim(temp)[2])
{
    #plot(density(as.numeric(as.character(temp[,i]))))
    sm.density.compare(as.numeric(temp[,i]), tgeClass.f, xlab=colnames(temp)[i])
    #title(main="Human Ovarian Cancer")
    colfill<-c(2:(2+length(levels(tgeClass.f))))
    #legend(locator(1), levels(tgeClass.f), fill=colfill)
   
}

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
temp$TGE_Class[grep("_",temp$TGE_Class)] <- "isoform"
temp_m = melt( temp[,2:dim(temp)[2]] , id="TGE_Class")
features = as.character(unique(temp_m$variable))

plot1 = list();
for(i in 1:length(features))
{
        temp_m_t = temp_m[ temp_m[,"variable"] == features[i], ]
        dp <- ggplot(temp_m_t, aes(value, fill = TGE_Class ) ) + geom_density(alpha = 0.5) + labs(title = features[i] ) + theme_bw()
        plot1[[i]] = dp;
}
library(grid)
library(gridExtra)
multiplot(plotlist = plot1, cols = 3)

pdf("C:/Users/shyama/Dropbox/PITDB/distribution.pdf")
invisible(lapply(plot1, print))
dev.off()
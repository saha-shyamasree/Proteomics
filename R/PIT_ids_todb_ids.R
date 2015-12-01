###this code reads PIT identifcation results and the blast result.
###Replaces the PIT Ids with the Uniprot/database ids..

source("D:/Code/Proteomics/R/RLib.R")
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

proteinGroupFiltered<- function(proteinGrp,rev,peptide,pepThreshold,revStr,conStr)
{
    proteinGrp=as.matrix(proteinGrp)
    print(dim(proteinGrp))
    if(rev==1)
    {
        proteinGrp=proteinGrp[-(grep(revStr,proteinGrp[,'protein.accession'])),]
        print(dim(proteinGrp))
        proteinGrp=proteinGrp[-(grep(conStr,proteinGrp[,'protein.accession'])),]
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
identifiedBlastResult<-function(Mat,B,out)
{
    ##B should hold blast result in CSV format
    ##Mat is protein identification result from MSGF+
    BIdentified=B[which(B[,'query_name'] %in% Mat[,'description']),]
    write.csv(BIdentified,file=out, quote = FALSE,row.names = FALSE)
}
readMats<- function(f1,f2,d1,d2,rev=1,peptide=1,pepThreshold=1,upper=0.1,revStr,contStr)
{
    #print(f1)
    #print(f2)
    #print(d1)
    #print(d2)
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold,revStr,contStr)
    blastDB=blastFilterEval(blastFilterMatch(blast(f2,d2),1),upper)
    print("Mat1")
    print(dim(Mat1))
    print("BLAST")
    print(dim(blastDB))
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    
    write.csv(Mat1,file=paste(d1,f1,"DBIds.csv",sep=""))
    identifiedBlastResult(Mat1, blastDB, paste(d2,"Identified/",f2,sep=""))
}

readMat<- function(f1,f2,d1,d2,rev=1,peptide=1,pepThreshold=1,upper=0.1,revStr,contStr)
{
    #print(f1)
    #print(f2)
    #print(d1)
    #print(d2)
    Mat1=proteinGroupFiltered(proteinGroup(f1,d1),rev,peptide,pepThreshold,revStr,contStr)
    blastDB=blastFilterEval(blastFilterMatch(blast(f2,d2),1),upper)
    print("Mat1")
    print(dim(Mat1))
    print("BLAST")
    print(dim(blastDB))
    Mat1[,'protein.accession']=replaceIds(Mat1,blastDB)
    list('Mat1'=Mat1,'blast'=blastDB)
}

sampleWiseProteinList<-function(dir1,dir2,rev=1,peptide=1,pepThreshold=1,upper=0.1,revStr,contStr)
{
    dirs=list.dirs(path=dir1,full.names = FALSE)
    protMat=NULL
    mat=list()
    blastDBs=list()
    known=list()
    kVar=list()
    iso=list()
    isovar=list()
    dirs=dirs[dirs!=""]
    dirs=dirs[dirs!="G30"]
    coln=c()
    for(i in dirs)
    {
        coln=append(coln,paste(i,c("Peptide","SAP","ISO","ISOVAR"),sep="_"))
    }

    for(i in 1:length(dirs))
    {
        print(i)
        d=dirs[i]
        file1=paste(d,".assemblies.fasta.transdecoder.pep+fdr+th+grouping+prt.csv",sep="")
        file2=paste(d,".assemblies.fasta.transdecoder.pep.xml2.csv",sep="")
        print(file1)
        d=paste(dir1,"/",d,"/",sep="")
        #print(d)
        #print(file2)
        #print(dir2)
        tempMat=readMat(file1,file2,d,dir2,rev,peptide,pepThreshold,upper,revStr,contStr)
        mat[[i]]=tempMat$Mat1
        blastDBs[[i]]=tempMat$blastDB
        baseName=gsub("\\.csv$","",file2)
        #known[[i]]=read.csv(paste(dir2,baseName,"_known.csv",sep=""),header=T)
        kVar[[i]]=read.csv(paste(dir2,baseName,"_knownVar.csv",sep=""),header=T)
        iso[[i]]=read.csv(paste(dir2,baseName,"_iso.csv",sep=""),header=T)
        isovar[[i]]=read.csv(paste(dir2,baseName,"_isoVar.csv",sep=""),header=T)
    }
    
    uniq_vals = c()
    for(i in 1:length(mat))
    {
        uniq_vals = append(uniq_vals,mat[[i]][,'protein.accession'])
    }
    uniq_vals = unique(uniq_vals)
    uniq_vals=gsub(" OS=.*","", gsub("sp\\||tr\\|","", uniq_vals))
    uniq_uniprot_vals=uniq_vals[-grep("^(asmbl_)",uniq_vals)]
    uniq_uniprot_vals=matrix(unlist(strsplit(uniq_uniprot_vals,c("\\|"))), ncol = 2, byrow = TRUE)
    
    ##Matrisome proteins
    matrisome=read.table("D:/Doc/Oliver/matrisom.txt",sep="\t",header=TRUE)
    matrisomeProt=unlist(strsplit(as.character(matrisome$UniProt_IDs),c("\\:")))
    canonicalUniprot=uniq_uniprot_vals
    canonicalUniprot[,1]=gsub("-\\d+","",canonicalUniprot[,1])
    idx=which( canonicalUniprot[,1] %in% matrisomeProt)
    uniq_uniprot_vals[idx,1]
    
    wrtMat = matrix( "", length(idx), length(dirs)*4 )
    rownames(wrtMat) = uniq_uniprot_vals[idx,1]
    colnames(wrtMat) = coln
    
    wrtList = list()
    
    set.seed(1)##find out whats the use of it.
    
    for(i in 1:length(mat))
    {
        uniprotMat=as.data.frame(mat[[i]][-grep("^(asmbl_)", mat[[i]][,'protein.accession']),])
        uniprotMat[,'protein.accession'] = gsub(" OS=.*","", gsub("sp\\||tr\\|","", uniprotMat[,'protein.accession']))
        idMat=matrix(unlist(strsplit(uniprotMat[,'protein.accession'],c("\\|"))), ncol = 2, byrow = TRUE)
        colnames(idMat)=c('UniprotID','protein_description')
        uniprotMatExt=cbind(uniprotMat,idMat)
        type=c("SAP","ISO","ISOVAR")
        uniprotMatExt[,type]='N'
        uniprotMatExt[which(uniprotMatExt[,'description'] %in% kVar[[i]][,'ORF.Id']),'SAP']='Y'
        uniprotMatExt[which(uniprotMatExt[,'description'] %in% iso[[i]][,'ORF.Id']),'ISO']='Y'
        uniprotMatExt[which(uniprotMatExt[,'description'] %in% isovar[[i]][,'ORF.Id']),'ISOVAR']='Y'
        matIdx=which(uniprotMatExt[,'UniprotID'] %in% uniq_uniprot_vals[idx,1])
        vals_mat=uniprotMatExt[matIdx,c('UniprotID','distinct.peptide.sequences','SAP','ISO','ISOVAR')]
        wrtList[[i]]=vals_mat 
        #vals_df=data.frame(UniprotID=vals_mat[,1],peptide=as.numeric(vals_mat[,2]),SAP=vals_mat[,'SAP'],ISO=vals_mat[,'ISO'],ISOVAR=vals_mat[,'ISOVAR'])
        #vals_matr=as.matrix(aggregate(.~UniprotID,data=vals_mat,FUN=paste, collapse = ";"))
        #vals_matr= vals_mat %>% group_by(UniprotID) %>% summarise_each(funs(toString)) 
        
        u_vals = as.character(unique(vals_mat[,1]))
        for(k in 1:length(u_vals))
        {
            temp = vals_mat[ vals_mat[,1] == u_vals[k], ]
            temp_vals = c()
            for( l in 2:dim(temp)[2] )
            {
                temp_vals [l-1] = paste(temp[,l],collapse=";")
            }
            temp_vals = append( u_vals[k] , temp_vals )
            if( k == 1 )
            {
                vals_matr =  temp_vals       
            }
            else{
                vals_matr = rbind( vals_matr,temp_vals)
            }
        }
        rownames(vals_matr)=vals_matr[,1]
        wrtMat[vals_matr[,1],((i - 1) * 4 + 1):((i * 1) * 4)] =  vals_matr[,2:5]
    }
    statMat=matrix(0,3,8)
    colnames(statMat)=dirs
    for(i in 1:length(dirs))
    {
        d=dirs[i]
        print(d)
        print(paste("dimension:",dim(wrtList[[i]])))
        print(paste("peptide Count:",sum(as.numeric(wrtList[[i]][,2]))))
        print(paste("SAP:",length(which(wrtList[[i]][,3]=="Y"))))
        statMat[1,i]=length(which(wrtList[[i]][,3]=="Y"))
        print(paste("ISO:",length(which(wrtList[[i]][,4]=="Y"))))
        statMat[2,i]=length(which(wrtList[[i]][,4]=="Y"))
        print(paste("ISOVAR:",length(which(wrtList[[i]][,5]=="Y"))))
        statMat[3,i]=length(which(wrtList[[i]][,5]=="Y"))
        print(paste("ECM Protein:",length(unique(wrtList[[i]][,1]))))
    }
    
    write.table(wrtMat,paste(dir1,"/wrtMat.tsv", sep=""), sep = "\t" , row.names=T, quote=FALSE)

}
ddply(vals_df,"UniprotID",numcolwise(paste))

df = data.frame(a=c("X","Y","X","Y","Z"),b=c(1,2,2,4,5),c=c(10,12,13,1,4))
aggregate(b+c ~ a, data=df, FUN=mean)

rev=1
peptide=1
pepThreshold=1
upper=0.000000000000000000000000000001
revStr="XXX_"
contStr="CONT_"
dir1="D:/data/Oliver/PASA"
dir2="D:/data/Oliver/Blast/Identified/Location/"

compMat=wrtMat
> grep("*_Peptide",colnames(compMat))
[1]  1  5  9 13 17 21 25 29
> colN=colnames(compMat)
> heatMat=compMat[,colN[grep("*_Peptide",colN)]]
> head(heatMat)
mode(heatMat)="numeric"
colnames(heatMat)=gsub("_Peptide","",colnames(heatMat))

colnames(statMatMelt)=c("Type","Sample","Count")
ggplot(statMatMelt, aes(x = Sample, y = Count, fill = Type)) + 
 scale_fill_manual(values = c("gray20","grey50","grey70")) +
 xlab("Sample") + ylab("ORF count")+ ggtitle("PIT identified ORFs similar to the ECM proteins") +
 geom_bar(stat = "identity")
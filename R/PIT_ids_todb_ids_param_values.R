###These are the parameters for PIT_ids_todb_ids.R

rev=1
peptide=1
pepThreshold=1
upper=0.000000000000000000000000000001
revStr="XXX_"
contStr="CONT_"
dir1="D:/data/Oliver/PASA"
dir2="D:/data/Oliver/Blast/"
dirs=list.dirs(path=dir1,full.names = FALSE)
for(d in dirs)
{
    if(d!="" && d!="G30")
    {
        file1=paste(d,".assemblies.fasta.transdecoder.pep+fdr+th+grouping+prt.csv",sep="")
        file2=paste(d,".assemblies.fasta.transdecoder.pep.xml2.csv",sep="")
        print(file1)
        d=paste(dir1,"/",d,"/",sep="")
        print(d)
        print(file2)
        print(dir2)
        readMats(file1,file2,d,dir2,rev,peptide,pepThreshold,upper,revStr,contStr)
    }
}

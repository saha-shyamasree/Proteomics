#!/bin/sh
#$ -cwd		#Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=240:0:0    # Request 24 hour runtime
#$ -l h_vmem=12G      # Request 32GB RAM

#source ~/.profile

#rna=$1
#Trinity=$2
#pasa=$1 #'/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'
#sample='G10'
#Trinity=$2 #'/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/'
#echo $Trinity

##Claude Rachel Samples

rdir='/data/SBCS-BessantLab/shyama/Data/Rachel_Claude/RNA/'
trinity='/data/SBCS-BessantLab/shyama/Data/Rachel_Claude/Trinity/'
for dir in "$rdir"*
do
    if [ -d "$dir" ];then
        subDir=$(basename "$dir")
        #echo $subDir
        
        if [ ! -d "$trinity$subDir" ];then
            mkdir "$trinity$subDir"
        fi
#<<COMMENT
        for file in "$dir"/*_1.fastq.gz
        do
            if [ -f $file ]; then
                #echo $file
                file2=$(echo $file | sed -e "s/_1.fastq.gz/_2.fastq.gz/g")
                #echo $file2
                sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
                echo $sample
                echo "Trinity --seqType fq --JM 24G --left $file --right $file2 --SS_lib_type RF --CPU 4 --trimmomatic --normalize_reads --output $trinity$subDir"/"$sample --full_cleanup" | qsub -cwd -V -l h_vmem=40G -l h_rt=180:0:0 -q smserial.q -pe smp 4
                #break
            fi
        done
#COMMENT
    fi
    #break
done

<<COMMENT
## remove the trinity directories
dir='/data/SBCS-BessantLab/shyama/Data/Claude/RNA/'
trinity='/data/SBCS-BessantLab/shyama/Data/Claude/Trinity/'
for file in "$dir"*_1.fastq.gz
do
    if [ -f $file ]; then
        if [[ $file != *"253"* ]] && [[ $file != *"254"* ]] && [[ $file != *"255"* ]]
        then
            sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
            echo $sample
            echo "rm -R $trinity$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
        else
            echo $file
        fi
    fi
done
COMMENT

<<COMMENT

dir='/data/SBCS-BessantLab/shyama/Data/Claude/RNA/'
trinity='/data/SBCS-BessantLab/shyama/Data/Claude/Trinity/'
for file in "$dir"*_1.fastq.gz
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1.fastq.gz/_2.fastq.gz/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
        echo $sample
        if [[ $file == *"253"* ]] || [[ $file == *"254"* ]] || [[ $file == *"255"* ]]
        then
            echo "Trinity --seqType fq --JM 24G --left $file --right $file2 --SS_lib_type RF --CPU 4 --trimmomatic --normalize_reads --output $trinity$sample --full_cleanup" | qsub -cwd -V -l h_vmem=40G -l h_rt=180:0:0 -q smserial.q -pe smp 4
        fi
    fi
done
COMMENT

<<COMMENT
## remove the trinity directories
dir='/data/SBCS-BessantLab/shyama/Data/Claude/RNA/'
trinity='/data/SBCS-BessantLab/shyama/Data/Claude/Trinity/'
for file in "$dir"*_1.fastq.gz
do
    if [ -f $file ]; then
        if [[ $file != *"253"* ]] && [[ $file != *"254"* ]] && [[ $file != *"255"* ]]
        then
            sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
            echo $sample
            echo "rm -R $trinity$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
        else
            echo $file
        fi
    fi
done
COMMENT
<<COMMENT
Trinity='/data/SBCS-BessantLab/shyama/Data/Claude/Trinity/'
DIR='/data/SBCS-BessantLab/shyama/Data/Claude/PASA/'
for file in "$Trinity"*.fasta
do
    sample=$(basename $file | sed -e "s/.Trinity.fasta//g")
    samp=$(echo $sample | sed -e 's/-/_/g')

    #if [ ! -f "$DIR$sample/$samp".assemblies.fasta ]; then
    #To run transdecoder
    if [ -f "$DIR$sample/$samp".assemblies.fasta ]; then
    	echo $DIR$samp/$samp.assemblies.fasta
        echo $samp
    	qsub pasa_pipeline_apocrita.sh -d $DIR -f $file -g /data/SBCS-BessantLab/shyama/Data/Fasta/hg38.fa -a /data/SBCS-BessantLab/shyama/Data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 -s 36
    #sh pasa_pipeline_apocrita.sh -d $DIR -f $file -g /data/SBCS-BessantLab/shyama/Data/Fasta/hg38.fa -a /data/SBCS-BessantLab/shyama/Data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 -c 1
    fi
    #break
done
COMMENT


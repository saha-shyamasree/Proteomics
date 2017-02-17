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

dir='/data/SBCS-BessantLab/shyama/Data/Prabhakar/RNA/'
<<COMMENT
for file in "$dir"*R1_001.fastq.gz
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/R1_001/R2_001/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_R1_001.fastq.gz//g")
        echo $sample
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --SS_lib_type RF --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Prabhakar/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done


## remove the trinity directories
dir='/data/SBCS-BessantLab/shyama/Data/Prabhakar/RNA/'
for file in "$dir"
do
    if [ -f $file ]; then
        sample=$(basename $file | sed -e "s/_R1_001.fastq.gz//g")
        echo "rm -R /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
    fi
done

COMMENT

Trinity='/data/SBCS-BessantLab/shyama/Data/Prabhakar/Trinity/'
DIR='/data/SBCS-BessantLab/shyama/Data/Prabhakar/PASA/'
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



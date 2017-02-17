#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=60:0:0    # Request 24 hour runtime
#$ -l h_vmem=25G      # Request 1GB RAM

fasta=$1
sample=$2
fastq1=$3
fastq2=$4
#out=$5

rsem-prepare-reference --bowtie2 -p 8 $fasta $sample

rsem-calculate-expression -p 16 --bowtie2 --paired-end $fastq1 $fastq2 $sample "$sample"
		

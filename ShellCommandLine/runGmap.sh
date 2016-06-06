#!/bin/sh
#$ -cwd         #Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=240:0:0    # Request 24 hour runtime
#$ -l h_vmem=32G 

genomeDir=$1
genomeDB=$2
inFile=$3
outFile=$4

gmap -D $genomeDir -d $genomeDB -t 8 -f 2 -n 0 -x 50 -B 5 --gff3-add-separator=0 -Y --min-identity=0.5 -O --intronlength=300000 $inFile > $outFile


#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=24:0:0    # Request 24 hour runtime
#$ -l h_vmem=1G      # Request 1GB RAM

makeblastdb -in /data/scratch/btw796/Data/HumanAdeno/fasta/Homo_sapiens.GRCh38.ncrna.fa -dbtype nucl -out /data/scratch/btw796/Data/HumanAdeno/fasta/Homo_sapiens.GRCh38.ncrna.aa
tblastn -db /data/scratch/btw796/Data/HumanAdeno/fasta/Homo_sapiens.GRCh38.ncrna.aa -query /data/scratch/btw796/Data/HumanAdeno/fasta/trinityPITORF_no_blast_match.fasta -gapopen 6 -gapextend 2 -out /data/scratch/btw796/Data/HumanAdeno/fasta/trinityPITORF_no_blast_out.xml -outfmt 5

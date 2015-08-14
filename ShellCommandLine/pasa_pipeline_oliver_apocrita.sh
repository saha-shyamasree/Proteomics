#!/bin/sh
#$ -cwd		#Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=240:0:0    # Request 24 hour runtime
#$ -l h_vmem=120G      # Request 32GB RAM

#source ~/.profile
DIR='/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'
sample='G10'
Trinity='/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/'

## make a new directory with sample name
#if [ ! -d "$DIR""$sample" ]; then
#	mkdir $DIR$sample
#fi
#cp $Trinity$sample'.Trinity.fasta' $DIR$sample
#cp ~/Prog/PASApipeline-master/human_adeno_data/alignAssembly.config $DIR$sample
#sed -i -- s/human_adeno_mydb_pasa/$sample/ $DIR$sample/alignAssembly.config
cd $DIR$sample
#seqclean $DIR$sample'/'$sample'.Trinity.fasta'
#/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/hg38.fa -t $sample'.Trinity.fasta.clean' -T -u $sample'.Trinity.fasta' --ALIGNERS blat,gmap --MAX_INTRON_LENGTH 300000 --ALT_SPLICE --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 16
## after it failed due to path to hg38, removing -C so that it does not try to creat databases
/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s 9 -R -g /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/hg38.fa -t $sample'.Trinity.fasta.clean' -T -u $sample'.Trinity.fasta' --ALIGNERS blat,gmap --MAX_INTRON_LENGTH 300000 --ALT_SPLICE --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 16

#../misc_utilities/pasa_gff3_validator.pl Homo_sapiens.GRCh38.78_Protein_Coding_Chr.gff3

#../scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g hg38.fa -P Homo_sapiens.GRCh38.78_Protein_Coding_Chr.gff3

#../scripts/Launch_PASA_pipeline.pl -c annotCompare.config -A -g hg38.fa -t  human_adenovirus_trinity_assembled_transcripts.fasta.clean --CPU 30

#blat hg38.fa blat_out_dir/partition.0.fa -q=rna -dots=100  -maxIntron=300000 -out=pslx -ooc=11.ooc blat_out_dir/partition.0.fa.pslx

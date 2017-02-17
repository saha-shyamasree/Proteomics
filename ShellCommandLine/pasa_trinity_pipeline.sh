#!/bin/sh
#$ -cwd		#Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=240:0:0    # Request 24 hour runtime
#$ -l h_vmem=12G      # Request 32GB RAM

#source ~/.profile
<<COMMENT
rna=$1
Trinity=$2
pasa=$1 #'/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'
#sample='G10'
#Trinity=$2 #'/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/'
echo $Trinity

dir='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
for file in "$dir"*.fastq*
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
        echo $sample
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done

## remove the trinity directories
dir='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
for file in "$dir"G7*_1.* "$dir"G8*_1* "$dir"G9*_1*
do
    if [ -f $file ]; then
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
        echo "rm -R /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
    fi
done


for file in "$Trinity"*.fasta*
do
    if [ -f $file ]; then
        #echo $file
        sample=$(basename $file | sed -e "s/.fasta.*//g") 
	#if [ "$sample" = 'G42' ] || [ "$sample" = 'G54' ] || [ "$sample" = 'G82' ] || [ "$sample" = 'G124' ]; then
	#if [ "$sample" = 'G69' ]; then #step 25
	#if [ "$sample" = 'G58' ] || [ "$sample" = 'G72' ] || [ "$sample" = 'G85' ]; then #step 24
	#if [ "$sample" = 'G93' ] || [ "$sample" = 'G94' ]; then #step 23
	#if [ "$sample" = 'G30' ]; then
		echo $sample
	
		#if [ ! -d "$DIR""$sample" ]; then
        	#	mkdir $DIR$sample
		#fi
		#cp $Trinity$sample'.Trinity.fasta' $DIR$sample
		#cp ~/Prog/PASApipeline-master/human_adeno_data/alignAssembly.config $DIR$sample
		#sed -i -- s/human_adeno_mydb_pasa/$sample/ $DIR$sample/alignAssembly.config
		#cd $DIR$sample
		#echo "seqclean $DIR$sample/$sample.Trinity.fasta" | qsub -cwd -V -l h_vmem=8G -l h_rt=72:0:0
		#echo "/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/hg38.fa -t $sample.Trinity.fasta.clean -T -u $sample.Trinity.fasta --ALIGNERS blat,gmap --MAX_INTRON_LENGTH 300000 --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 8" | qsub -cwd -V -l h_vmem=60G -l h_rt=140:0:0
		#echo "/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s 23 -R -g /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/hg38.fa -t $sample.Trinity.fasta.clean -T -u $sample.Trinity.fasta --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 8" | qsub -cwd -V -l h_vmem=100G -l h_rt=100:0:0
		#echo "/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -g /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/hg38.fa -t $sample.Trinity.fasta.clean -T -u $sample.Trinity.fasta --ALT_SPLICE --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 4" | qsub -cwd -V -l h_vmem=16G -l h_rt=48:0:0
	#fi
    fi
done
COMMENT

## remove the trinity directories
dir='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
trinity='/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/'
for file in "$dir"*_1.fastq.gz
do
    if [ -f $file ]; then
        #if [[ $file != *"253"* ]] && [[ $file != *"254"* ]] && [[ $file != *"255"* ]]
        #then
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
        #echo $sample
        if [ -d $trinity$sample ]; then
            echo "rm -R $trinity$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
            echo $sample
        else
            echo $file
        fi
    fi
done
## make a new directory with sample name
#if [ ! -d "$DIR""$sample" ]; then
#	mkdir $DIR$sample
#fi
#cp $Trinity$sample'.Trinity.fasta' $DIR$sample
#cp ~/Prog/PASApipeline-master/human_adeno_data/alignAssembly.config $DIR$sample
#sed -i -- s/human_adeno_mydb_pasa/$sample/ $DIR$sample/alignAssembly.config
#cd $DIR$sample
#seqclean $DIR$sample'/'$sample'.Trinity.fasta'
#/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/hg38.fa -t $sample'.Trinity.fasta.clean' -T -u $sample'.Trinity.fasta' --ALIGNERS blat,gmap --MAX_INTRON_LENGTH 300000 --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 8
## after it failed due to path to hg38, removing -C so that it does not try to creat databases
#/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s 29 -R -g /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/hg38.fa -t $sample'.Trinity.fasta.clean' -T -u $sample'.Trinity.fasta' --ALIGNERS blat,gmap --MAX_INTRON_LENGTH 300000 --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 8

#/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl --ALT_SPLICE -c alignAssembly.config -g /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/hg38.fa -t $sample'.Trinity.fasta.clean' --CPU 8

#../misc_utilities/pasa_gff3_validator.pl Homo_sapiens.GRCh38.78_Protein_Coding_Chr.gff3

#../scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g hg38.fa -P Homo_sapiens.GRCh38.78_Protein_Coding_Chr.gff3

#../scripts/Launch_PASA_pipeline.pl -c annotCompare.config -A -g hg38.fa -t  human_adenovirus_trinity_assembled_transcripts.fasta.clean --CPU 30

#blat hg38.fa blat_out_dir/partition.0.fa -q=rna -dots=100  -maxIntron=300000 -out=pslx -ooc=11.ooc blat_out_dir/partition.0.fa.pslx

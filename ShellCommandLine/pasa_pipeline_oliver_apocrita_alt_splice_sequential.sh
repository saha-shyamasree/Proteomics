#!/bin/sh
#$ -cwd		#Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=240:0:0    # Request 24 hour runtime
#$ -l h_vmem=120G      # Request 32GB RAM

#source ~/.profile

DIR='/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'
#sample='G10'
Trinity='/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/'
for file in "$Trinity"G*.Trinity.fasta
do
    if [ -f $file ]; then
        #echo $file
        sample=$(basename $file | sed -e "s/.Trinity.fasta//g") 
	if [ "$sample" = 'G11' ] || [ "$sample" = 'G15' ] || [ "$sample" = 'G17' ] || [ "$sample" = 'G29a' ] || [ "$sample" = 'G75' ] || [ "$sample" = 'G102' ] || [ "$sample" = 'G138' ]; then
	#[ "$sample" = 'G10' ] ||

		echo $sample
		cd $DIR$sample
		/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -g /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/hg38.fa -t $sample.Trinity.fasta.clean -T -u $sample.Trinity.fasta --ALT_SPLICE --TRANSDECODER --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 /data/home/btw796/Prog/PASApipeline-master/human_adeno_data/Homo_sapiens.GRCh38.78_Protein_Coding.gff3 --CPU 8
	fi
    fi
done



DIR='/data/SBCS-BessantLab/shyama/Data/Oliver/Tophat/'
#sample='G10'
rna='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
index='/data/SBCS-BessantLab/shyama/Data/Fasta/hg38'
for file in "$rna"G*_1.fastq.gz
do
    if [ -f $file ]; then
        #echo $file
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
        #if [ "$sample" = 'G10' ] || [ "$sample" = 'G102' ] || [ "$sample" = 'G45' ] || [ "$sample" = 'G138' ] || [ "$sample" = 'G15' ] || [ "$sample" = 'G54' ] || [ "$sample" = 'G36' ] || [ "$sample" = 'G92' ]; then
        if [ "$sample" = 'G29a' ] || [ "$sample" = 'G30' ] || [ "$sample" = 'G43' ] || [ "$sample" = 'G51' ] || [ "$sample" = 'G79' ]; then
	#if [ "$sample" = 'G33a' ]; then
                echo $sample
		folder=$DIR$sample
                #if [ ! -d "$folder" ]; then
		       
                #       mkdir $folder
                #fi
                #cp $Trinity$sample'.Trinity.fasta' $DIR$sample
                #cp ~/Prog/PASApipeline-master/human_adeno_data/alignAssembly.config $DIR$sample
                #sed -i -- s/human_adeno_mydb_pasa/$sample/ $DIR$sample/alignAssembly.config
                cd $folder
		#left="$rna$sample"_1.fastq.gz
		#right="$rna$sample"_2.fastq.gz
		#echo "tophat2 -o $folder --library-type fr-firststrand --b2-very-fast -p 8 $index $left,$right" | qsub -cwd -V -l h_vmem=30G -l h_rt=100:0:0
		#echo "cufflinks -g /data/SBCS-BessantLab/shyama/Data/GTF/Homo_sapiens.GRCh38.78.gtf -p 8 --library-type fr-firststrand -o $folder $folder/accepted_hits.bam" | qsub -cwd -V -l h_vmem=16G -l h_rt=100:0:0
		echo "cuffdiff "
        fi
    fi
done

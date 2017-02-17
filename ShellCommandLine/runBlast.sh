
#echo "makeblastdb -in /data/SBCS-BessantLab/shyama/Data/Oliver/fasta/uniprot-proteomeUP000005640.fasta -dbtype prot -out /data/SBCS-BessantLab/shyama/Data/Oliver/fasta/uniprot-proteomeUP000005640.aa" | qsub -cwd -V -l h_vmem=20G -l h_rt=100:0:0

DIR='/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'

Trinity='/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/'
#xmlDIR='/data/SBCS-BessantLab/shyama/Data/Oliver/BLAST/xml/'
xmlDIR='/data/SBCS-BessantLab/shyama/Data/Oliver/BLAST/identified/xml/'
out='/data/SBCS-BessantLab/shyama/Data/Oliver/BLAST/identified/csv/'
for file in "$Trinity"G*.Trinity.fasta
do
    if [ -f $file ]; then
        sample=$(basename $file | sed -e "s/.Trinity.fasta//g")
        #if [ "$sample" = 'G8' ] || [ "$sample" = 'G10' ] || [ "$sample" = 'G11' ] || [ "$sample" = 'G15' ] || [ "$sample" = 'G17' ] || [ "$sample" = 'G29a' ] || [ "$sample" = 'G51' ] || [ "$sample" = 'G75' ] || [ "$sample" = 'G102' ] || [ "$sample" = 'G138' ]; then
#        if [ "$sample" = 'G30' ]; then
		if [ -f $DIR$sample/$sample.assemblies.fasta.transdecoder.pep.identified.fasta ]; then
		    echo $sample
<<COMMENT
			cd $DIR$sample
			if [ -f $sample.assemblies.fasta.transdecoder.pep.identified.nostart.fasta ]; then
				echo "exist"
			else
				cp $sample.assemblies.fasta.transdecoder.pep.identified.fasta $sample.assemblies.fasta.transdecoder.pep.identified.nostar.fasta
				sed -i 's/\*//g' $sample.assemblies.fasta.transdecoder.pep.identified.nostar.fasta
			fi

			echo "blastp -db /data/SBCS-BessantLab/shyama/Data/Oliver/fasta/uniprot-proteomeUP000005640.aa -query $sample.assemblies.fasta.transdecoder.pep.identified.nostar.fasta -gapopen 6 -gapextend 2 -out $xmlDIR$sample.assemblies.fasta.transdecoder.pep.identified.xml -outfmt 5 -evalue 1 -matrix BLOSUM80" | qsub -cwd -V -l h_vmem=8G -l h_rt=50:0:0

			cd /data/home/btw796/Code2/SYBARIS/Python
			python3.4 contigStat.py $xmlDIR$sample.assemblies.fasta.transdecoder.pep.identified.xml $out
COMMENT
			#python UniProteinLocation.py /data/SBCS-BessantLab/shyama/Data/
		fi
    fi
done


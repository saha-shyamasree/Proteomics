##download file from cluster
##Mouse
#<<COMMENT
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/'
blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/blast/'
annotOut='/mnt/c/Users/shyama/Dropbox/PITDB/MouseNelsonBay/Annotation/'
RSEMC='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/RSEM/'
RSEMD='/mnt/g/Bristol/Mouse/RSEM/'
for sampleDir in $pasa*
do
	if [ -d $sampleDir ]; then
		sample=$(basename $sampleDir)
        echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$RSEMC$sample/$sample.isoforms.results.identified.tsv $RSEMD$sample.isoforms.results.identified.tsv"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$blast$sample/$sample.assemblies.fasta.transdecoder.pep_annotation.csv $annotOut"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 ./"
 		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$blast$sample/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf ./"
 		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified_peptide.gff3 ./"
	fi
done
#COMMENT
##Bat Nelson Bay
#<<COMMENT2
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/'
blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/blast/'
annotOut='/mnt/c/Users/shyama/Dropbox/PITDB/BatNelsonBay/Annotation/'
RSEMC='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/RSEM/'
RSEMD='/mnt/g/Bristol/Bat/NelsonBay/RSEM/'
for sampleDir in $pasa*
do
	if [ -d $sampleDir ]; then
		sample=$(basename $sampleDir)
        echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$RSEMC$sample/$sample.isoforms.results.identified.tsv $RSEMD$sample.isoforms.results.identified.tsv"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$blast$sample/$sample.assemblies.fasta.transdecoder.pep_annotation.csv $annotOut"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 ./"
 		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$blast$sample/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf ./"
 		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified_peptide.gff3 ./"
	fi
done

##Mosquito
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/'
blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/blast/'
annotOut='/mnt/c/Users/shyama/Dropbox/PITDB/Mosquito/Annotation/'
sample='aedes'
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$RSEMC/$sample.isoforms.results.identified.tsv $RSEMD"
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$blast/$sample.assemblies.fasta.transdecoder.pep_annotation.csv $annotOut"
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 ./"
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$blast/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf ./"
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified_peptide.gff3 ./"

##Human Adeno
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/'
blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/'
annot='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/'
sample='human_adeno_mydb_pasa'
annotOut='/mnt/c/Users/shyama/Dropbox/PITDB/HumanAdeno/Annotation/'
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$annot/$sample.assemblies.fasta.transdecoder.pep_annotation.csv $annotOut"
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 ./"
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$blast/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf ./"
#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified_peptide.gff3 ./"

#COMMENT2

##Oliver's data
pasa='/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'
annot='/data/SBCS-BessantLab/shyama/Data/Oliver/BLAST/identified/csv/'
ident='/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/'
annotOut='/mnt/c/Users/shyama/Dropbox/PITDB/OliverData/Annotation/'
RSEMC='/data/SBCS-BessantLab/shyama/Data/Oliver/RSEM/'
RSEMD='/mnt/g/Oliver/RSEM/'
for sampleDir in $annot*
do
	if [ -d $sampleDir ]; then
		sample=$(basename $sampleDir)
        echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$RSEMC$sample/$sample.isoforms.results.identified.tsv $RSEMD$sample.isoforms.results.identified.tsv"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$annot$sample/$sample.assemblies.fasta.transdecoder.pep_annotation.csv $annotOut"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 ./"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.identified.fasta ./"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.pep.identified.fasta ./"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$ident$sample/$sample+fdr+th+grouping+prt_filtered.csv ./"
		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$ident$sample/$sample+fdr+th+grouping_filtered.csv ./"
 		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$annot$sample/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf ./"
 		#echo "scp btw796@frontend1.apocrita.hpc.qmul.ac.uk:$pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified_peptide.gff3 ./"
	fi
done


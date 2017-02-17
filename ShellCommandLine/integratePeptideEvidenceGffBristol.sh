###This code calls UniProteinLocation, IdentifyProteinIsoformSAP, peptideEvidence for Bristol data (Mouse, Bat nelson bay)
#source activate quantdb3.4

python='/data/home/btw796/Code2/Proteomics/Python'
integrateCode=$python/integratePeptideEvidenceInGFF3.py

###Mouse Nelson Bay
#<<COMMENT
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/identification/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/'
for sampleDir in $pasa*
do
    if [ -d $sampleDir ]; then
        sample=$(basename $sampleDir)
    	echo $sample
    	if [ -f $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv ]; then
    		#if [ "$sample" == "7_C5LH8ACXX_CAGATC_L005" ] || [ "$sample" == "8_C5LH8ACXX_ACTTGA_L005" ] || [ "$sample" == "9_C5LH8ACXX_GATCAG_L005" ] || [ "$sample" == "10_C5LH8ACXX_TAGCTT_L005" ]; then
			echo "python $integrateCode --psm $ident$sample/$sample+fdr+th+grouping_filtered.csv --gff3 $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 --proteins $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv" | qsub -cwd -l h_vmem=32G -l h_rt=24:0:0
			#fi
		fi
	fi
done


###Bat Nelson Bay

ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/'
for sampleDir in $pasa*
do
    if [ -d $sampleDir ]; then
        sample=$(basename $sampleDir)
        echo $sample
        #if [ "$sample" == "10_C5LH8ACXX_TAGCTT_L005" ] || [ "$sample" == "14_C5LH8ACXX_AGTTCC_L005" ] || [ "$sample" == "17_C5LH8ACXX_GTCCGC_L005" ]; then
		echo "python $integrateCode --psm $ident$sample/$sample+fdr+th+grouping_filtered.csv --gff3 $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 --proteins $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv" | qsub -cwd -l h_vmem=64G -l h_rt=24:0:0
		#fi
	fi
done

ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/identification/PASA/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/'
sample='aedes'

echo $sample	
echo "python $integrateCode --psm $ident$sample+fdr+th+grouping_filtered.csv --gff3 $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 --proteins $ident$sample+fdr+th+grouping+prt_filtered.csv" | qsub -cwd -l h_vmem=64G -l h_rt=24:0:0
#COMMENT
<<COMMENT1
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/'
sample='human_adeno_mydb_pasa'

echo $sample	
echo "python $integrateCode --psm $ident$sample+fdr+th+grouping_filtered.csv --gff3 $pasa$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 --proteins $ident$sample+fdr+th+grouping+prt_filtered.csv" | qsub -cwd -l h_vmem=64G -l h_rt=24:0:0
COMMENT


###This code calls UniProteinLocation, IdentifyProteinIsoformSAP, peptideEvidence for Bristol data (Mouse, Bat nelson bay)
#source activate quantdb3.4

python='/data/home/btw796/Code2/Proteomics/Python'
integrateCode=$python/integratePeptideEvidenceInGFF3.py

###Human 
ident='/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/'
pasa='/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'
for sampleDir in $pasa*
do
    if [ -d $sampleDir ]; then
        sample=$(basename $sampleDir)
    	echo $sample
    	if [ -f $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv ]; then
    		#if [ "$sample" == "G102" ]; then
				echo "python $integrateCode --psm $ident$sample/$sample+fdr+th+grouping_filtered.csv --gff3 $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3_identified.gff3 --proteins $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv" | qsub -cwd -l h_vmem=32G -l h_rt=24:0:0
			#fi
		fi
	fi
done




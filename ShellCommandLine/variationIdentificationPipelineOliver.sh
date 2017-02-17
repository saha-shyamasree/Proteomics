###This code calls UniProteinLocation, IdentifyProteinIsoformSAP, peptideEvidence for Bristol data (Mouse, Bat nelson bay)
#source activate quantdb3.4

python='/data/home/btw796/Code2/Proteomics/Python/IsoformsSAP'
uniloc=$python/UniProteinLocation.py
isoSAP=$python/IdentifyProteinIsoformSAP.py
pepEvd=$python/peptideEvidence.py
annotCode=$python/annotationMatrix.py
preliAnnot=$python/PreliminaryProteinAnnotationForPITDBV2.py
###Mouse Nelson Bay

blast='/data/SBCS-BessantLab/shyama/Data/Oliver/BLAST/identified/csv/'
uniLoc='/data/SBCS-BessantLab/shyama/Data/uniprot-Homo+sapiens960613_04_2015.tsv'
out='/data/SBCS-BessantLab/shyama/Data/Oliver/BLAST/identified/csv/'
ident='/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/'
pasa='/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/'

for sampleDir in $blast*
do
    if [ -d $sampleDir ]; then
        sample=$(basename $sampleDir)
        echo $sample
    fi
    if [ -f $sampleDir ]; then
    	sample=$(basename $sampleDir | sed -e "s/.assemblies.fasta.transdecoder.pep.identified.xml2.csv//g")
    	echo $sample
    	if [ ! -d $blast$sample ]; then
    		mkdir $blast$sample
    		mv $sampleDir $blast$sample
    	fi
    fi    
		#python3.4 $uniloc -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.xml2.csv -u $uniLoc -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
	
		#python3.4 $isoSAP -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv -k $sampleDir/$sample.assemblies.fasta.transdecoder.pep_known.csv -s $sampleDir/$sample.assemblies.fasta.transdecoder.pep_knownVar.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -i $sampleDir/$sample.assemblies.fasta.transdecoder.pep_iso.csv -j $sampleDir/$sample.assemblies.fasta.transdecoder.pep_isoVar.csv > $python/$sample.log
		#echo "python3.4 $pepEvd -p $ident$sample/$sample+fdr+th+grouping_filtered.csv -s $ident$sample/$sample+fdr+th+grouping_filtered_variation.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf"
		#echo "python3.4 $pepEvd -p $ident$sample/$sample+fdr+th+grouping_filtered.csv -s $ident$sample/$sample+fdr+th+grouping_filtered_variation.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf" | qsub -cwd -l h_vmem=16G -l h_rt=120:0:0
		#python3.4 $annotCode -p $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv -g $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3
		python3.4 $preliAnnot -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
done


###This code calls UniProteinLocation, IdentifyProteinIsoformSAP, peptideEvidence for Bristol data (Mouse, Bat nelson bay)
#source activate quantdb3.4

python='/data/home/btw796/Code2/Proteomics/Python/IsoformsSAP'
uniloc=$python/UniProteinLocation.py
isoSAP=$python/IdentifyProteinIsoformSAP.py
pepEvd=$python/peptideEvidence.py
annotCode=$python/annotationMatrix.py
preliAnnot=$python/PreliminaryProteinAnnotationForPITDBV2.py
###Mouse Nelson Bay

blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/blast/'
uniLoc='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/fasta/uniprot-taxonomy10090+NelsonBay.tsv'
out='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/blast/'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/identification/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/'
for sampleDir in $blast*
do
    if [ -d $sampleDir ]; then
        sample=$(basename $sampleDir)
        echo $sample
        
		#python3.4 $uniloc -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.xml2.csv -u $uniLoc -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
	
		#python3.4 $isoSAP -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv -k $sampleDir/$sample.assemblies.fasta.transdecoder.pep_known.csv -s $sampleDir/$sample.assemblies.fasta.transdecoder.pep_knownVar.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -i $sampleDir/$sample.assemblies.fasta.transdecoder.pep_iso.csv -j $sampleDir/$sample.assemblies.fasta.transdecoder.pep_isoVar.csv > $python/$sample.log
		#echo "python3.4 $pepEvd -p $ident$sample/$sample+fdr+th+grouping_filtered.csv -s $ident$sample/$sample+fdr+th+grouping_filtered_variation.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf"
		#echo "python3.4 $pepEvd -p $ident$sample/$sample+fdr+th+grouping_filtered.csv -s $ident$sample/$sample+fdr+th+grouping_filtered_variation.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf" | qsub -cwd -l h_vmem=16G -l h_rt=120:0:0
		#python3.4 $annotCode -p $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv -g $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3
		python3.4 $preliAnnot -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
	fi
done


###Bat Nelson Bay

blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/blast/'
uniLoc='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/fasta/uniprot-organismPteropus+alecto9402+NelsonBay.tsv'
out='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/blast/'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/'
for sampleDir in $blast*
do
    if [ -d $sampleDir ]; then
        sample=$(basename $sampleDir)
        #echo $sample
		#python3.4 $uniloc -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.xml2.csv -u $uniLoc -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
		
		#python3.4 $isoSAP -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv -k $sampleDir/$sample.assemblies.fasta.transdecoder.pep_known.csv -s $sampleDir/$sample.assemblies.fasta.transdecoder.pep_knownVar.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -i $sampleDir/$sample.assemblies.fasta.transdecoder.pep_iso.csv -j $sampleDir/$sample.assemblies.fasta.transdecoder.pep_isoVar.csv > $python/$sample.log
		
		#echo "python3.4 $pepEvd -p $ident$sample/$sample+fdr+th+grouping_filtered.csv -s $ident$sample/$sample+fdr+th+grouping_filtered_variation.csv -v $sampleDir/$sample.assemblies.fasta.transdecoder.pep.vcf -o $sampleDir/$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf" | qsub -cwd -l h_vmem=16G -l h_rt=120:0:0
		#python3.4 $annotCode -p $ident$sample/$sample+fdr+th+grouping+prt_filtered.csv -g $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3
		#break
		python3.4 $preliAnnot -b $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
	fi
done


blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/blast/'
uniLoc='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/fasta/uniprot-proteomeUP000008820.tsv'
out='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/blast/'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/identification/PASA/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/'
sample='aedes'
#python3.4 $uniloc -b $blast$sample.assemblies.fasta.transdecoder.pep.identified.xml2.csv -u $uniLoc -o $blast$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
		
#python3.4 $isoSAP -b $blast$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv -k $blast$sample.assemblies.fasta.transdecoder.pep_known.csv -s $blast$sample.assemblies.fasta.transdecoder.pep_knownVar.csv -v $blast$sample.assemblies.fasta.transdecoder.pep.vcf -i $blast$sample.assemblies.fasta.transdecoder.pep_iso.csv -j $blast$sample.assemblies.fasta.transdecoder.pep_isoVar.csv > $python/$sample.log
		
#echo "python3.4 $pepEvd -p $ident$sample+fdr+th+grouping_filtered.csv -s $ident$sample+fdr+th+grouping_filtered_variation.csv -v $blast$sample.assemblies.fasta.transdecoder.pep.vcf -o $blast$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf" | qsub -cwd -l h_vmem=16G -l h_rt=120:0:0
#python3.4 $annotCode -p $ident/$sample+fdr+th+grouping+prt_filtered.csv -g $pasa$sample/$sample.assemblies.fasta.transdecoder.genome.gff3
python3.4 $preliAnnot -b $blast$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv


blast='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/'
uniLoc='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/uniprot-Homo+sapiens960613_04_2015.tsv'
out='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/blast/'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/'
sample='human_adeno_mydb_pasa'
#python3.4 $uniloc -b $blast$sample.assemblies.fasta.transdecoder.pep.identified.xml2.csv -u $uniLoc -o $blast$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv
		
#python3.4 $isoSAP -b $blast$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv -k $blast$sample.assemblies.fasta.transdecoder.pep_known.csv -s $blast$sample.assemblies.fasta.transdecoder.pep_knownVar.csv -v $blast$sample.assemblies.fasta.transdecoder.pep.vcf -i $blast$sample.assemblies.fasta.transdecoder.pep_iso.csv -j $blast$sample.assemblies.fasta.transdecoder.pep_isoVar.csv > $python/$sample.log
		
#echo "python3.4 $pepEvd -p $ident$sample+fdr+th+grouping_filtered.csv -s $ident$sample+fdr+th+grouping_filtered_variation.csv -v $blast$sample.assemblies.fasta.transdecoder.pep.vcf -o $blast$sample.assemblies.fasta.transdecoder.pep_pepEvd.vcf" | qsub -cwd -l h_vmem=16G -l h_rt=120:0:0
#python3.4 $annotCode -p $ident$sample+fdr+th+grouping+prt_filtered.csv -g $pasa$sample.assemblies.fasta.transdecoder.genome.gff3
python3.4 $preliAnnot -b $blast$sample.assemblies.fasta.transdecoder.pep.identified.loc.csv


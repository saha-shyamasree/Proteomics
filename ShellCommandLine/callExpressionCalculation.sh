##This script call expression calculation code.

usage()
{
	echo "usage: $0 [-d baseDir] "
}

while getopts d: opt
do
    case $opt in
      (d)  baseDir="$OPTARG";;
      (\?)
      	  usage
	  exit;;
    esac
done
shift `expr $OPTIND - 1`

#echo $baseDir

if [ -z "$baseDir" ]; then
	usage
	exit
fi

rna=$baseDir"RNA/"
count=0

##uncooment this line when running peptideEvidenceIsoformStats.py
#echo -e "Sample\tTotal Isoforms\tSingle-sided\tDouble-sided\tShortened\tPeptideEvd(TGEObs)\tPeptideEvd(Variation)\tUniqueEvd(TGEObs)\tUniqueEvd(Variation)\tJuntionPeptide"
##for Prab's data file path is like PC3M-Nsi1_S1_R1_001.fastq.gz, i.e. paired ends are R1 and R2
for fastq1 in "$rna"*R1_001.fastq* #"$rna"*1.fastq*
do
	echo "$fastq1"
    ##Bristol and Oliver's data
	#sample=$(basename $fastq1| sed -r "s/_.?1\.fastq(.*)?//g")
	#fastq2=`echo "$fastq1" | sed -r "s/1\.fastq/2.fastq/g"`
    ##For Prab's data
    sample=$(basename $fastq1| sed -r "s/_R1_001\.fastq\.gz//g")
    fastq2=`echo "$fastq1" | sed -r "s/_R1_001\.fastq\.gz/_R2_001.fastq.gz/g"`
	#echo "$fastq2"
	#echo $sample
    ##Prabhakar's sample name contains '-', that pasa did not like, hence it has been replaced by '_'
    pSample=`echo "$sample" | sed -r "s/-/_/g"`
    echo $pSample
    
    if [ -e "$baseDir"RSEM/"$sample"/"$pSample".assemblies.fasta.trinity.fasta ]; then
        echo $pSample
        qsub expressionCalculation.sh "$baseDir"RSEM/"$sample"/"$pSample".assemblies.fasta.trinity.fasta "$baseDir"RSEM/"$sample"/"$pSample" "$fastq1" "$fastq2"
    fi
    #if [ -e "$baseDir"RSEM/"$sample"/"$sample".transcript.bam ] && [ -e "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv ]; then
        #echo $sample
        #cd ../Python/IsoformsSAP/
        #cd /data/home/btw796/Code2/Proteomics/Python
        #python FPKM_of_identified_ORFs.py -f "$baseDir"RSEM/"$sample"/"$sample".isoforms.results -a "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv -o "$baseDir"RSEM/"$sample"/"$sample".isoforms.results.identified.tsv
    #fi
    ##following code is to run peptideEvidenceIsoform.py . We have a different line for Oliver's data because the file is different. in Future that should not happen.
    #if [ -e "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf ]; then
        #cd ../Python/IsoformsSAP/
        #echo "python peptideEvidenceIsoforms.py -b "$baseDir"blast/"$sample"/"$sample".assemblies.fasta.transdecoder.pep.identified.loc.csv -a "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv -p "$baseDir"PITDB/PSMs-Peptides-ORFs/"$sample"+fdr+th+grouping_filtered.csv" | qsub -cwd -l h_rt=10:0:0 -l h_vmem=16G
        #echo "File exists"
    #else
        #if [ -e "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv ]; then
            #cd ../Python/IsoformsSAP/
            #echo "python peptideEvidenceIsoforms.py -b "$baseDir"blast/"$sample"/"$sample".assemblies.fasta.transdecoder.pep.identified.loc.csv -a "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv -p "$baseDir"PITDB/PSMs-Peptides-ORFs/"$sample"+fdr+th+grouping_filtered.csv" | qsub -cwd -l h_rt=10:0:0 -l h_vmem=8G
            ##Oliver's data
            #echo "python peptideEvidenceIsoforms.py -b "$baseDir"BLAST/identified/csv/"$sample"/"$sample".assemblies.fasta.transdecoder.pep.identified.loc.csv -a "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv -p "$baseDir"PITDB/PSMs-Peptides-ORFs/"$sample"+fdr+th+grouping_filtered.csv" | qsub -cwd -l h_rt=10:0:0 -l h_vmem=16G
        #fi
    #fi
    #if [ -e "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf ]; then
        #echo "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf
    #    python peptideEvidenceIsoformsStats.py -v "$baseDir"PITDB/Annotation/"$sample".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf
    #fi
done


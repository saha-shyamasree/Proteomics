
genomeDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data"
genomeDB="hg38.fa.gmap"
inFile="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G10/G10.assemblies.fasta"
outFile="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G10/G10.assemblies.gmap_gene.gff3"
qsub runGmap.sh $genomeDir $genomeDB $inFile $outFile

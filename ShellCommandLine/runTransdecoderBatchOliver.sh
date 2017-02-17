##This code run runTransdecoder.sh in batch mode.

PASA="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/"

for path in "$PASA"*
do
	#echo $path
	sample=$(basename $path)
	#echo "$PASA$sample"/"$sample".assemblies.fasta
	if [ -f "$path"/"$sample".assemblies.fasta ] && [ ! -f "$path"/"$sample".assemblies.fasta.transdecoder.pep ]; then
		echo $sample
		sh runTransdecoder.sh $path $sample
	fi
done

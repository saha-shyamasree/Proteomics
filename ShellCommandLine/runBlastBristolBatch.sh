#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=120:0:0    # Request 240 hour runtime
#$ -l h_vmem=16G      # Request 4GB RAM
###This code create BLAST database from Uniprot proteome fasta file and BLAST the identified ORFs (passed as argument) agaist it.

usage()
{
	echo "usage: $0 [-u proteomefasta] [-p orffolder] [-o outfolder]"
}


while getopts u:p:o: opt
do
    case $opt in
      (u)  uni="$OPTARG";;
      (p)  DIR="$OPTARG";;
      (o)  out="$OPTARG";;
      (\?)
      	  usage
	  exit;;
    esac
done
shift `expr $OPTIND - 1`

if [ -z "$uni" ]; then
	usage
	exit
fi

if [ -z "$DIR" ]; then
        usage
        exit
fi

if [ -z "$out" ]; then
        usage
        exit
fi

echo $uni
echo $DIR
echo $out

proteome=$(basename $uni | sed -e "s/.fasta//g") 
#echo $proteome

fastaPath=$(dirname $uni)
#echo $fastaPath

if [ -f $fastaPath/$proteome.aa.pin ]; then
	echo "BLAST DB exist"
else
	echo "makeblastdb -in $uni -dbtype prot -out $fastaPath/$proteome.aa"
	#makeblastdb -in $uni -dbtype prot -out $fastaPath/$proteome.aa
fi

for sampleDir in $DIR*
do
    if [ -d $sampleDir ]; then
        sample=$(basename $sampleDir)
        #echo $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.fasta
        #echo "\n\n"
		if [ -f $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.fasta ]; then
			echo $sample
			echo "qsub runBlastBristol.sh -u $uni -p $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.fasta -o $out/$sample"
			qsub runBlastBristol.sh -u $uni -p $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.fasta -o $out/$sample
			#sh runBlastBristol.sh -u $uni -p $sampleDir/$sample.assemblies.fasta.transdecoder.pep.identified.fasta -o $out/$sample
		fi
    fi
done

##HUman adeno stdin.o3758143


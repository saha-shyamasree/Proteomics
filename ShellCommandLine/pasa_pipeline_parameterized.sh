#!/bin/sh
#$ -cwd         #Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=120:0:0    # Request 24 hour runtime
#$ -l h_vmem=120G      # Request 32GB RAM
usage()
{
	echo "usage: $0 [-d outdirectory*] [-f trinityfasta*] [-g genomefasta*] [-a annotation] *options are mandatory"
}

unset DIR
unset file
unset genome
unset annotation
while getopts d:f:g:a: opt
do
    case $opt in
      (d)  DIR="$OPTARG";;
      (f)  file="$OPTARG";;
      (g)  genome="$OPTARG";;
      (a)  annotation="$OPTARG";;
      (\?)
      	  usage
	  exit;;
    esac
done
shift `expr $OPTIND - 1`
#echo $file
if [ -z "$file" ]; then
	usage
	exit
fi

if [ -z "$DIR" ]; then
        usage
        exit
fi

if [ -z "$genome" ]; then
        usage
        exit
fi


if [ -f $file ]; then
	echo $file
        sample=$(basename $file | sed -e "s/.fasta//g") 
	echo $sample
	fname=$(basename $file)
	if [ ! -d "$DIR""$sample" ]; then
        	mkdir $DIR$sample
        fi
        cp $file $DIR$sample
        cp /data/SBCS-BessantLab/shyama/Data/alignAssembly.config $DIR$sample
        sed -i -- s/human_adeno_mydb_pasa/$sample/ $DIR$sample/alignAssembly.config

	cd $DIR$sample
	echo $DIR$sample/$fname
	seqclean $DIR$sample/$fname
	if [ -z "$annotation" ]; then
		/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $genome -t $DIR$sample/$fname.clean -T -u $DIR$sample/$fname --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --stringent_alignment_overlap 30.0 --gene_overlap 50.0 --CPU 16
	else
		/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $genome -t $DIR$sample/$fname.clean -T -u $DIR$sample/$fname --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --stringent_alignment_overlap 30.0 --gene_overlap 50.0 -L --annots_gff3 $annotation --CPU 16
	fi
	/data/home/btw796/Prog/PASApipeline-master/scripts/pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta $sample.assemblies.fasta --pasa_transcripts_gff3 $sample.pasa_assemblies.gff3	
fi


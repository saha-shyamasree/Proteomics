#!/bin/sh
#$ -cwd         #Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=120:0:0    # Request 24 hour runtime
#$ -l h_vmem=120G      # Request 32GB RAM
usage()
{
	echo "usage: $0 [-d outdirectory] [-f trinityfasta] [-g genomefasta] [-a annotation in GFF3 format] [-c database create flag [0/1(default)]] [-s if a job failed after creating the database, provide this option to run the pipeline from given index.]"
}

unset DIR
unset file

while getopts d:f:g:a:c:s: opt
do
    case $opt in
      (d)  DIR="$OPTARG";;
      (f)  file="$OPTARG";;
      (g)  genome="$OPTARG";;
      (a)  annotation="$OPTARG";;
      (c)  create="$OPTARG";;
      (s)  step="$OPTARG";;
      (\?)
      	  usage
	  exit;;
    esac
done
shift `expr $OPTIND - 1`

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

if [ -z "$annotation" ]; then
        usage
        exit
fi

if [ -z "$create" ]; then
	if [ -z "$step" ]; then
        	usage
        	exit
        fi
fi


echo $file
echo $genome
echo $annotation
echo $create
echo $step


if [ -f $file ]; then
	echo $file
    sample=$(basename $file | sed -e "s/.Trinity.fasta//g") 
	echo $sample
	if [ ! -d "$DIR""$sample" ]; then
        	mkdir $DIR$sample
        fi
        cd $DIR$sample
        def=1
        samp=$(echo $sample | sed -e 's/-/_/g')
		echo $samp
        if [ -n "$create" ]; then
		if [ "$create" -eq "$def" ]
		then
			echo "create data base"
			cp $file $DIR$sample
			cp /data/SBCS-BessantLab/shyama/Data/alignAssembly.config $DIR$sample
			
			sed -i -- s/human_adeno_mydb_pasa/$samp/ $DIR$sample/alignAssembly.config		
			seqclean $DIR$sample/$sample'.Trinity.fasta'
			/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $genome -t $DIR$sample/$sample'.Trinity.fasta.clean' -T -u $DIR$sample/$sample'.Trinity.fasta' --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --stringent_alignment_overlap 30.0 --gene_overlap 50.0 --CPU 16
			#/data/home/btw796/Prog/PASApipeline-master2/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g $genome -t $DIR$sample/$sample'.Trinity.fasta.clean' -T -u $DIR$sample/$sample'.Trinity.fasta' --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --stringent_alignment_overlap 30.0 --gene_overlap 50.0 --CPU 16
			##run transdecoder
			/data/home/btw796/Prog/PASApipeline-master/scripts/pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta $sample.assemblies.fasta --pasa_transcripts_gff3 $sample.pasa_assemblies.gff3	
			#/data/home/btw796/Prog/PASApipeline-master2/PASApipeline-master/scripts/pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta $sample.assemblies.fasta --pasa_transcripts_gff3 $sample.pasa_assemblies.gff3
		fi
        fi
        
        if [ -n "$step" ]; then
        	echo "Step:"$step
        	#/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s $step -R -g $genome -t $DIR$sample/$sample'.Trinity.fasta.clean' -T -u $DIR$sample/$sample'.Trinity.fasta' --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --stringent_alignment_overlap 30.0 --gene_overlap 50.0 --CPU 8
        	#/data/home/btw796/Prog/PASApipeline-master2/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s $step -R -g $genome -t $DIR$sample/$sample'.Trinity.fasta.clean' -T -u $DIR$sample/$sample'.Trinity.fasta' --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --stringent_alignment_overlap 30.0 --gene_overlap 50.0 --CPU 8
        	##run transdecoder
			/data/home/btw796/Prog/PASApipeline-master/scripts/pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta $samp.assemblies.fasta --pasa_transcripts_gff3 $samp.pasa_assemblies.gff3	
			#/data/home/btw796/Prog/PASApipeline-master2/PASApipeline-master/scripts/pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta $sample.assemblies.fasta --pasa_transcripts_gff3 $sample.pasa_assemblies.gff3
        fi
        
	##Only to run alternative splicing
	#/data/home/btw796/Prog/PASApipeline-master/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -g $genome -t $DIR$sample/$sample'.Trinity.fasta.clean' -T -u $DIR$sample/$sample'.Trinity.fasta' --ALIGNERS blat,gmap --ALT_SPLICE --MAX_INTRON_LENGTH 300000 --stringent_alignment_overlap 30.0 --gene_overlap 50.0 --CPU 8
		
fi


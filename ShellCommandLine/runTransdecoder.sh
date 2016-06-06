DIR=$1
sample=$2
#sample should be the DB name used during PASA run.
echo $DIR
echo $sample

cd $DIR
echo "/data/home/btw796/Prog/PASApipeline-master/scripts/pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta $sample.assemblies.fasta --pasa_transcripts_gff3 $sample.pasa_assemblies.gff3" | qsub -cwd -V -l h_vmem=8G -l h_rt=24:0:0


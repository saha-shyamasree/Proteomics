
inD=$1
#outD=$2

for file in "$inD"*.fasta
do
	echo $file
	echo "makeblastdb -in $file -dbtype prot -out $file.aa" | qsub -cwd -V -l h_vmem=10G -l h_rt=100:0:0
done

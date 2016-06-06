
inF=$1
outF=$2
#cd /data/home/btw796/Code/Proteomics/Python
cd /data/home/btw796/Code2/SYBARIS/Python
#python3.4 contigStat.py $inF $outF
for file in "$inF"*.xml
do
	echo $file
	echo "python3.4 contigStat.py $file $outF" | qsub -cwd -l h_vmem=4G -l h_rt=50:0:0
done

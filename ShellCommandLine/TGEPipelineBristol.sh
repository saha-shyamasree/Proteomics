##### This code assembles David Mathew's Mouse and Bat infected with nelson bay virus data. #####
##Mouse
dir='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/RNA/'
for file in "$dir"*_R1.fastq
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_R1.fastq//g") 
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done


##Bat, nelson bay
dir='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/RNA/'
for file in "$dir"*_R1.fastq
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_R1.fastq//g") 
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done


for file in "$dir"12*_R1.fastq "$dir"14*_R1.fastq "$dir"15*_R1.fastq "$dir"17*_R1.fastq
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_R1.fastq//g") 
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done

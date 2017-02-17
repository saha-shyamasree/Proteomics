# This is the pipeline to sample wise assembles RNA-Seq data using Trinity. Before carrying out the assembly, for paired end stranded data
# one sample is mapped to the genome, if available, to find out library type. Mapped reads need to be manually observed to identify the read
# orientation.
############################################################################################################################################
###################################################### This is for model organisms #########################################################
############################################################################################################################################
## 76, 81, 88
## Get the genome from Ensembl
<<COMMENT
wget ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz /data/SBCS-BessantLab/shyama/Data/Oliver/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
## Build bowtie2 index from the refernec genome
echo "bowtie2-build /data/SBCS-BessantLab/shyama/Data/Oliver/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa /data/SBCS-BessantLab/shyama/Data/Oliver/fasta/Homo_sapiens.GRCh38.dna.toplevel" | qsub -cwd -V -l h_vmem=20G -l h_rt=48:0:0
## Map first 3000 reads to the genome for strand specific library type identification
echo "bowtie2 -q -u 3000 --no-mixed --qc-filter -x /data/SBCS-BessantLab/shyama/Data/Oliver/fasta/Homo_sapiens.GRCh38.dna.toplevel -1 /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G10_1.fastq.gz -2 /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G10_2.fastq.gz -S /data/SBCS-BessantLab/shyama/Data/Oliver/libraryTypeExp/G10.sam" | qsub -cwd -V -l h_vmem=16G -l h_rt=24:0:0

## After this step I load the sam file in IGV. The sam file is sorted and index using igvtools from IGV and viewed to identify the strand specific library type, i.e.
## forward-reverse or reverse-forward.

## Following command should take each sample and assemble them.
## Run trinity to do both de novo and genome guided assembly
## de novo
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G10_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G10_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G10 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G11_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G11_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G11 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G15_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G15_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G15 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G17_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G17_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G17 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G29a_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G29a_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G29a --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G30_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G30_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G30 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G33a_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G33a_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G33a --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G36_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G36_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G36 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G42_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G42_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G42 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G43_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G43_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G43 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G45_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G45_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G45 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G51_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G51_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G51 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G54_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G54_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G54 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G57_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G57_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G57 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G58_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G58_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G58 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G67_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G67_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G67 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
echo "Trinity --seqType fq --JM 30G --left /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G69_1.fastq.gz --right /data/SBCS-BessantLab/shyama/Data/Oliver/RNA/G69_2.fastq.gz --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/G69 --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0

dir='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
for file in "$dir"G7*_1.* "$dir"G8*_1* "$dir"G9*_1*
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g") 
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done

dir='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
for file in "$dir"G1[0-9][0-9]_1.*
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g") 
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done
dir='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
for file in "$dir"G8[18]_1.* "$dir"G76_1.*
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g")
        echo $sample
        echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --SS_lib_type FR --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
    fi
done

## remove the trinity directories
dir='/data/SBCS-BessantLab/shyama/Data/Oliver/RNA/'
for file in "$dir"G7*_1.* "$dir"G8*_1* "$dir"G9*_1*
do
    if [ -f $file ]; then
        sample=$(basename $file | sed -e "s/_1.fastq.gz//g") 
        echo "rm -R /data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
    fi
done
COMMENT

DIR='/data/SBCS-BessantLab/shyama/Data/Oliver/PITDB/PSMs-Peptides-ORFs/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/'

code='/data/home/btw796/Code2/Proteomics/Python/standardSearchResultProcessing.py'
#cd code
for file in "$Trinity"*.Trinity.fasta
do
    sample=$(basename $file | sed -e "s/.Trinity.fasta//g")
    cp /data/SBCS-BessantLab/shyama/Data/Oliver/Identification/$sample/$sample.standard+fdr+th+grouping+prt.csv $DIR/$sample.standard+fdr+th+grouping+prt.csv
    cp /data/SBCS-BessantLab/shyama/Data/Oliver/Identification/$sample/$sample.standard+fdr+th+grouping.csv $DIR/$sample.standard+fdr+th+grouping.csv
    #mv $DIR$sample+fdr+th+grouping_filtered.csv $DIR$sample+fdr+th+grouping_filtered_bug.csv
	#python $code --protein $DIR$sample+fdr+th+grouping+prt_filtered.csv --peptide $DIR/$sample+fdr+th+grouping_filtered_bug.csv --protOut $DIR$sample+fdr+th+grouping+prt_filtered.csv --pepOut $DIR$sample+fdr+th+grouping_filtered.csv
    python $code --protein $DIR$sample.standard+fdr+th+grouping+prt.csv --peptide $DIR/$sample.standard+fdr+th+grouping.csv --protOut $DIR$sample.standard+fdr+th+grouping+prt_filtered.csv --pepOut $DIR$sample.standard+fdr+th+grouping_filtered.csv
	
done
##### This code assembles David Mathew's Mouse and Bat infected with nelson bay virus data. After that it runs 'runProteinIdentificationAndPostProcessing_cluster.py' to generate the protein sidentification commands#####
###It runs for Mosquito data as well.
###It runs Human Hendra virus as well.
##Mouse

dir='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/RNA/'
trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/Trinity/'
pasa='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/identification/'
<<COMMENT
for file in "$dir"*_R1.fastq
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_R1.fastq//g") 
        #echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
        #if [ -d $trinity$sample ]; then
        #    echo "rm -R $trinity$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
        #fi
        ##Merge ORFs of all samples to check whether the merged result is correct.
        if [ $sample != "1_C5LH8ACXX_ATCACG_L005" ]; then
            echo $sample
            cat $pasa$sample/$sample.assemblies.fasta.transdecoder.pep >> $pasa"allSample.pep"
        fi
    fi
done
COMMENT
##This code is to run MSGF on all the samples ORFs
sample="allSampleMouse"
python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf -o $ident$sample/$sample.mzid -d $pasa"allSampleMouse.pep" -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications.txt > run_$sample.cmds.sh
qsub run_$sample.cmds.sh
<<COMMENT

##Bat, nelson bay
dir='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/RNA/'
trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/'
for file in "$dir"*_R1.fastq
do
    if [ -f $file ]; then
        #echo $file
        file2=$(echo $file | sed -e "s/_1/_2/g")
        #echo $file2
        sample=$(basename $file | sed -e "s/_R1.fastq//g") 
        #echo "Trinity --seqType fq --JM 30G --left $file --right $file2 --CPU 6 --trimmomatic --normalize_reads --output /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/$sample --full_cleanup" | qsub -cwd -V -l h_vmem=80G -l h_rt=72:0:0
        echo "rm -R $trinity$sample" | qsub -cwd -V -l h_vmem=2G -l h_rt=72:0:0
    fi
done
COMMENT

<<COMMENT1

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
COMMENT1
<<COMMENT11
##Mouse NelsonBay
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/identification/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/Trinity/'
#genome='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/fasta/GCF_000325575.1_ASM32557v1_genomic.fasta'
codePath='/data/home/btw796/Code2/Proteomics/Python'
cd $codePath
for file in "$Trinity"*.Trinity.fasta
do
        sample=$(basename $file | sed -e "s/.Trinity.fasta//g")
        one="1_"
        two="2_"
        three="3_"
        four="4_"
        five="5_"
        six="6_"
        seven="7_"
        eight="8_"
        nine="9_"
        
        #if test "${file#*$one}" != "$file" || test "${file#*$two}" != "$file" || test "${file#*$three}" != "$file"
        #then
        #	if [ ! -d "$ident""$sample" ]; then
        #		mkdir $ident$sample
        #	fi
        #       echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_light.txt > $sample.cmds.txt"
        #        python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_light.txt > run_$sample.cmds.sh
        #fi
        
        #if test "${file#*$four}" != "$file" || test "${file#*$five}" != "$file" || test "${file#*$six}" != "$file"
        #then
        #	if [ ! -d "$ident""$sample" ]; then
        #		mkdir $ident$sample
        #	fi
        #	echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_medium.txt > $sample.cmds.txt"
        #        python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_medium.txt > run_$sample.cmds.sh
        #fi
        
        if test "${file#*$seven}" != "$file" || test "${file#*$eight}" != "$file" || test "${file#*$nine}" != "$file"
        then
        	if [ ! -d "$ident""$sample" ]; then
        		mkdir $ident$sample
        	fi
        	echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_heavy.txt > $sample.cmds.txt"
                python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_heavy.txt > run_$sample.cmds.sh

        fi
done

COMMENT11

##Bat NelsonBay 
<<COMMENT2
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/'
genome='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/fasta/GCF_000325575.1_ASM32557v1_genomic.fasta'
code='/data/home/btw796/Code2/Proteomics/ShellCommandLine/pasa_pipeline_apocrita.sh'
#cd code
for file in "$Trinity"*.Trinity.fasta
do
	if [ "$file" != '/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/10_C5LH8ACXX_TAGCTT_L005.Trinity.fasta' ]; then
		echo "sh $code -d $DIR -f $file -g $genome" | qsub -cwd -V -l h_vmem=120G -l h_rt=120:0:0
		#echo "sh $code -d $DIR -f $file -g $genome"
	fi
	#qsub -v d=$DIR,f=$Trinity,g=$genome $code
done
COMMENT2
<<COMMENT3
##Bat NelsonBay
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/'
genome='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/fasta/GCF_000325575.1_ASM32557v1_genomic.fasta'
codePath='/data/home/btw796/Code2/Proteomics/Python/'
cd $codePath

for file in "$Trinity"*.Trinity.fasta
do
	sample=$(basename $file | sed -e "s/.Trinity.fasta//g")
	ten="10_"
	eleven="11_"
	twelve="12_"
	thirteen="13_"
	fourteen="14_"
	fifteen="15_"
	sixteen="16_"
	seventeen="17_"
	eighteen="18_"
	if test "${file#*$ten}" != "$file" || test "${file#*$eleven}" != "$file" || test "${file#*$twelve}" != "$file"
	then
		if [ ! -d "$ident""$sample" ]; then
        		mkdir $ident$sample
        	fi
		echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/P.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_light.txt > run_$sample.cmds.sh"
		python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/P.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_light.txt > run_$sample.cmds.sh

        fi
        if test "${file#*$thirteen}" != "$file" || test "${file#*$fourteen}" != "$file" || test "${file#*$fifteen}" != "$file"
	then
		if [ ! -d "$ident""$sample" ]; then
        		mkdir $ident$sample
        	fi
		echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/P.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_medium.txt > run_$sample.cmds.sh"
		python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/P.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_medium.txt > run_$sample.cmds.sh

        fi
        if test "${file#*$sixteen}" != "$file" || test "${file#*$seventeen}" != "$file" || test "${file#*$eighteen}" != "$file"
	then
		if [ ! -d "$ident""$sample" ]; then
        		mkdir $ident$sample
        	fi
		echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/P.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_heavy.txt > run_$sample.cmds.sh"
		python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/P.mgf -o $ident$sample/$sample.mzid -d $DIR$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_SILAC_heavy.txt > run_$sample.cmds.sh

        fi
        #qsub -v d=$DIR,f=$Trinity,g=$genome $code
done
COMMENT3

<<COMMENT4
####Mosquito
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/aedes'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/identification/PASA'

sample='aedes'
codePath='/data/home/btw796/Code2/Proteomics/Python/'
cd $codePath

echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/MS/SliceAll.mgf -o $ident/$sample.mzid -d $DIR/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_mosquito.txt > run_$sample.cmds.sh"
python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/MS/SliceAll.mgf -o $ident/$sample.mzid -d $DIR/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_mosquito.txt > run_$sample.cmds.sh

COMMENT4

<<COMMENT5
## Bat Hendra
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/Hendra/PASA/Trinity_P_alecto'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/Hendra/identification/PASA'

sample='Trinity_P_alecto'
codePath='/data/home/btw796/Code2/Proteomics/Python/'
cd $codePath

echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/Hendra/MS/Slice.mgf -o $ident/$sample.mzid -d $DIR/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications.txt > run_$sample.cmds.sh"
python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/Hendra/MS/Slice.mgf -o $ident/$sample.mzid -d $DIR/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications.txt > run_$sample.cmds.sh



## Human Hendra

DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/hendra/PASA/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/hendra/Trinity/'
genome='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/fasta/GCF_000325575.1_ASM32557v1_genomic.fasta'
code='/data/home/btw796/Code2/Proteomics/ShellCommandLine/pasa_pipeline_apocrita.sh'

COMMENT5

####Human Adeno
<<COMMENT
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus'
ident='/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA'

sample='human_adeno_mydb_pasa.custom.All.standard'
codePath='/data/home/btw796/Code2/Proteomics/Python/'
cd $codePath

echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i /data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/MS/SliceAll.mgf -o $ident/$sample.mzid -d $DIR/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Bristol/modifications_mosquito.txt > run_$sample.cmds.sh"
python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i $DIR/MS/DM_from_raw.mgf -o $ident/$sample.mzid -d $DIR/PITDB/AminoAcids-or-ORFs-orTGEs/human_adeno.Custom.All.Standard.fasta -m $DIR/MS/modifications.txt > run_$sample.cmds.sh
qsub run_$sample.cmds.sh
COMMENT
<<COMMENT
##Bat NelsonBay,
#This is to only update the filtered peptide file that has already been produced, but was wrong due to the regex issue.
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PITDB/PSMs-Peptides-ORFs/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/'

code='/data/home/btw796/Code2/Proteomics/Python/standardSearchResultProcessing.py'
#cd code
for file in "$Trinity"*.Trinity.fasta
do
    sample=$(basename $file | sed -e "s/.Trinity.fasta//g")
    cp /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/$sample/$sample+fdr+th+grouping+prt_filtered.csv $DIR/$sample+fdr+th+grouping+prt_filtered.csv
    #mv $DIR/$sample+fdr+th+grouping_filtered.csv $DIR/$sample+fdr+th+grouping_filtered_bug.csv
	python $code --protein $DIR/$sample+fdr+th+grouping+prt_filtered.csv --peptide $DIR/$sample+fdr+th+grouping_filtered_bug.csv --protOut $DIR/$sample+fdr+th+grouping+prt_filtered.csv --pepOut $DIR/$sample+fdr+th+grouping_filtered.csv
	
done
COMMENT
<<COMMENT
## Mouse NelsonBay,
#This is to only update the filtered peptide file that has already been produced, but was wrong due to the regex issue.
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PITDB/PSMs-Peptides-ORFs/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/Trinity/'

code='/data/home/btw796/Code2/Proteomics/Python/standardSearchResultProcessing.py'
#cd code
for file in "$Trinity"*.Trinity.fasta
do
    sample=$(basename $file | sed -e "s/.Trinity.fasta//g")
    #cp /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/$sample/$sample+fdr+th+grouping+prt_filtered.csv $DIR/$sample+fdr+th+grouping+prt_filtered.csv
    mv $DIR$sample+fdr+th+grouping_filtered.csv $DIR$sample+fdr+th+grouping_filtered_bug.csv
	python $code --protein $DIR$sample+fdr+th+grouping+prt_filtered.csv --peptide $DIR/$sample+fdr+th+grouping_filtered_bug.csv --protOut $DIR$sample+fdr+th+grouping+prt_filtered.csv --pepOut $DIR$sample+fdr+th+grouping_filtered.csv
	
done
COMMENT
<<COMMENT
## Mosquito
#This is to only update the filtered peptide file that has already been produced, but was wrong due to the regex issue.
DIR='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PITDB/PSMs-Peptides-ORFs/'
Trinity='/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/Trinity/'

code='/data/home/btw796/Code2/Proteomics/Python/standardSearchResultProcessing.py'
#cd code
for file in "$Trinity"*Trinity.fasta
do
    sample=$(basename $file | sed -e "s/Trinity.fasta//g")
    #cp /data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/$sample/$sample+fdr+th+grouping+prt_filtered.csv $DIR/$sample+fdr+th+grouping+prt_filtered.csv
    mv $DIR$sample+fdr+th+grouping_filtered.csv $DIR$sample+fdr+th+grouping_filtered_bug.csv
	python $code --protein $DIR$sample+fdr+th+grouping+prt_filtered.csv --peptide $DIR/$sample+fdr+th+grouping_filtered_bug.csv --protOut $DIR$sample+fdr+th+grouping+prt_filtered.csv --pepOut $DIR$sample+fdr+th+grouping_filtered.csv
	
done
COMMENT
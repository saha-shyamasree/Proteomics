
inD=$1
#file=$inD"ORF.fasta"
for file in "$inD"*ORF.fasta
do
	echo "$file"
	for db in "$inD"*Jun2015.fasta
	do
		echo "$db"
		dbStr=""
		if [[ $db == *"all"* ]]
		then
			dbStr="All"
			echo "All"
			 echo "blastp -db $db.aa -query $file -gapopen 6 -gapextend 2 -out "$file"_"$dbStr".xml -outfmt 5 -evalue 1" | qsub -cwd -l h_vmem=4G -l h_rt=200:0:0
		fi
		if [[ $db == *"reviewed"* ]]
		then
			dbStr="reviewed"
			echo "Reviewed"
		fi	
		
		#echo "blastp -db $db.aa -query $file -gapopen 6 -gapextend 2 -out "$file"_"$dbStr".xml -outfmt 5 -evalue 1" | qsub -cwd -l h_vmem=4G -l h_rt=24:0:0 
	done	
done

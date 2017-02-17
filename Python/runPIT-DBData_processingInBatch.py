### This code is for running 'PIT-DBData_processing.py' in batch mode.
import os
import glob
import re

'''
orfDir="D:/data/Oliver/ORFs/"
transDir="D:/data/Oliver/transcripts/"
identified="D:/data/Oliver/PASA/"

samples=["G10","G11","G15","G17","G138","G75","G29a","G102"]

for s in samples:
    #command="python PIT-DBData_processing.py --ORFs "+orfDir+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+".assemblies.fasta --proteins "+identified+s+"/"+s+".assemblies.fasta.transdecoder.pep+fdr+th+grouping+prt.csvDBIds.csv --peptides "+identified+s+"/"+s+".assemblies.fasta.transdecoder.pep+fdr+th+grouping.csv --ORFsOutFile "+orfDir+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile "+transDir+s+".assemblies.fasta.identified.fasta --peptideOutFile "+identified+s+"/"+s+".assemblies.fasta.transdecoder.pep+fdr+th+grouping_filtered.csv"
    command=""
    #print(command)
    print("\n\n")
    #os.system(command)

##Mosquito

orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/aedes/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/aedes/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/identification/PASA/"

s="aedes"
command="python3.4 PIT-DBData_processing.py --ORFs "+orfDir+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+".assemblies.fasta --proteins "+identified+s+"+fdr+th+grouping+prt.csv --peptides "+identified+s+"+fdr+th+grouping.csv --ORFsOutFile "+orfDir+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile "+transDir+s+"assemblies.fasta.identified.fasta --proteinOutFile "+identified+s+"+fdr+th+grouping+prt_filtered.csv --peptideOutFile "+identified+s+"+fdr+th+grouping_filtered.csv"
#print(command)
print("\n\n")
#os.system(command)

##Bat Nelson Bay

trinityDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/"
orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/"

for f in glob.glob(trinityDir+"*.Trinity.fasta"):
	print(f)
	s=re.sub("\.Trinity\.fasta","",os.path.basename(f))
	print(s)
	command="python3.4 PIT-DBData_processing.py --ORFs "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+"/"+s+".assemblies.fasta --proteins "+identified+s+"/"+s+"+fdr+th+grouping+prt.csv --peptides "+identified+s+"/"+s+"+fdr+th+grouping.csv --ORFsOutFile "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile "+transDir+s+"/"+s+".assemblies.fasta.identified.fasta --proteinOutFile "+identified+s+"/"+s+"+fdr+th+grouping+prt_filtered.csv --peptideOutFile "+identified+s+"/"+s+"+fdr+th+grouping_filtered.csv"
	#print(command)
	print("\n")
	#os.system(command)


##Mouse nelson bay

trinityDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/Trinity/"
orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/identification/"

for f in glob.glob(trinityDir+"*.Trinity.fasta"):
	print(f)
	s=re.sub("\.Trinity\.fasta","",os.path.basename(f))
	print(s)
	command="python3.4 PIT-DBData_processing.py --ORFs "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+"/"+s+".assemblies.fasta --proteins "+identified+s+"/"+s+"+fdr+th+grouping+prt.csv --peptides "+identified+s+"/"+s+"+fdr+th+grouping.csv --ORFsOutFile "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile "+transDir+s+"/"+s+".assemblies.fasta.identified.fasta --proteinOutFile "+identified+s+"/"+s+"+fdr+th+grouping+prt_filtered.csv --peptideOutFile "+identified+s+"/"+s+"+fdr+th+grouping_filtered.csv"
	
	if os.path.isfile(identified+s+"/"+s+"+fdr+th+grouping+prt.csv"):
		#print(command)
		print("\n")
		#os.system(command)
'''	
##Human adeno-virus
'''
orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/"

s="human_adeno_mydb_pasa"
#command="python3.4 PIT-DBData_processing.py --ORFs "+orfDir+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+".assemblies.fasta --proteins "+identified+s+"+fdr+th+grouping+prt.csv --peptides "+identified+s+"+fdr+th+grouping.csv --ORFsOutFile "+orfDir+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile "+transDir+s+".assemblies.fasta.identified.fasta --proteinOutFile "+identified+s+"+fdr+th+grouping+prt_filtered.csv --peptideOutFile "+identified+s+"+fdr+th+grouping_filtered.csv"
#print(command)
print("\n\n")
#os.system(command)

commandGff="python IsoformsSAP/annotationMatrix.py -p "+identified+s+"+fdr+th+grouping+prt_filtered.csv -g "+orfDir+s+".assemblies.fasta.transdecoder.genome.gff3"
print(commandGff)
os.system(commandGff)
'''
##Oliver's data

trinityDir="/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/"
orfDir="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/"
transDir="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/"
identified="/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/"
processedSamples=['G10','G11','G15','G17','G29a','G75','G102','G138']

for f in glob.glob(trinityDir+"*.Trinity.fasta"):
	print(f)
	s=re.sub("\.Trinity\.fasta","",os.path.basename(f))
	print(s)
	command="python3.4 PIT-DBData_processing.py --ORFs "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+"/"+s+".assemblies.fasta --proteins "+identified+s+"/"+s+"+fdr+th+grouping+prt.csv --peptides "+identified+s+"/"+s+"+fdr+th+grouping.csv --ORFsOutFile "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile "+transDir+s+"/"+s+".assemblies.fasta.identified.fasta --proteinOutFile "+identified+s+"/"+s+"+fdr+th+grouping+prt_filtered.csv --peptideOutFile "+identified+s+"/"+s+"+fdr+th+grouping_filtered.csv"
	command2="python3.4 IsoformsSAP/annotationMatrix.py -p "+identified+s+"/"+s+"+fdr+th+grouping+prt.csv -g "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.genome.gff3"
	if os.path.isfile(identified+s+"/"+s+"+fdr+th+grouping+prt.csv"):
		if s not in processedSamples:
			#print(command)
			#os.system(command)
			print(command2)
			#print("\n")
			os.system(command2)


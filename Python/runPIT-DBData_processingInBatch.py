### This code is for running 'PIT-DBData_processing.py' in batch mode.
import os

orfDir="D:/data/Oliver/ORFs/"
transDir="D:/data/Oliver/transcripts/"
identified="D:/data/Oliver/PASA/"

samples=["G10","G11","G15","G17","G138","G75","G29a","G102"]

for s in samples:
    command="python PIT-DBData_processing.py --ORFs "+orfDir+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+".assemblies.fasta --proteins "+identified+s+"/"+s+".assemblies.fasta.transdecoder.pep+fdr+th+grouping+prt.csvDBIds.csv --peptides "+identified+s+"/"+s+".assemblies.fasta.transdecoder.pep+fdr+th+grouping.csv --ORFsOutFile "+orfDir+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsOutFile "+transDir+s+".assemblies.fasta.identified.fasta --peptideOutFile "+identified+s+"/"+s+".assemblies.fasta.transdecoder.pep+fdr+th+grouping_filtered.csv"
    print(command)
    print("\n\n")
    os.system(command)
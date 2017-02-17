##This code calls virusTrinscriptsExtractionAndORFPrediction.py
import os
import glob
import argparse

parser = argparse.ArgumentParser(description='This python code find failed transcripts for both the gmap and blat aligner. These are possible virus transcripts. This code then extracts the trinity trascript and save in a file for further processing.')
parser.add_argument("-b", "--base", nargs=1, required=True, help="full path of to the base directory of the dataset", metavar="PATH")
#parser.add_argument("-t", "--trinity", nargs=1, required=True, help="full path of to the trinity fasta", metavar="PATH")
#parser.add_argument("-p", "--pasa", nargs=1, required=True, help="full path of to the PASA fasta", metavar="PATH")
#parser.add_argument("-o", "--out", nargs=1, required=True, help="full path of to the output fasta file", metavar="PATH")

args = parser.parse_args()
print(args)
shellCode="/data/home/btw796/Code2/Proteomics/ShellCommandLine/expressionCalculation.sh"
pasaDir=args.base[0]+"PASA/"
for sample in os.listdir(pasaDir):
	#print("Sample"+sample)
	#if sample=="G69":
	if os.path.isdir(pasaDir+sample):
		print("Sample:"+sample)
		##Prabhakars sample name contains '-'. PASA had problem with that, hence that character has been replaced by '_'
		pSample=sample.replace("-","_")
		pasa=pasaDir+sample+"/"+pSample+".assemblies.fasta"
		bv=pasaDir+sample+"/"+pSample+".valid_blat_alignments.bed"
		bf=pasaDir+sample+"/"+pSample+".failed_blat_alignments.bed"
		gv=pasaDir+sample+"/"+pSample+".valid_gmap_alignments.bed"
		gf=pasaDir+sample+"/"+pSample+".failed_gmap_alignments.bed"
		trinity=pasaDir+sample+"/"+sample+".Trinity.fasta"
		rsem=args.base[0]+"RSEM/"
		if not os.path.isdir(rsem):
			os.makedirs(rsem)
		if not os.path.isdir(rsem+sample+"/"):
			os.makedirs(rsem+sample+"/")
		if os.path.exists(pasa):
			print("pasa exists")
			os.system("python virusTranscriptExtractionAndORFPrediction.py -l "+bv+" -m "+gv+" -b "+bf+" -g "+gf+" -t "+trinity+" -p "+pasa+" -o " \
		+rsem+sample+"/"+sample+"_missing.fasta")
		else:
			print(pasa+"does not exist")
		if os.path.exists(pasa+".trinity.fasta"):
			#print(pasa+".trinity.fasta exists" )
			os.system("mv "+pasa+".trinity.fasta "+rsem+sample+"/")
		#break

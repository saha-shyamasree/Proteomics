##this code extracts trinity transcript ids from SAMPLE_NAME.failed_blat_alignments.bed and SAMPLE_NAME.failed_gmap_alignments.bed.
##When a transcript id appears in both the files, it suggests that the transcript was not used for PASA assembly. Such transcripts
##are possible virus transcripts.

import os
import argparse
import pandas as pd
import glob
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def readFile(filename, sep, headerFlag):
	fileDFObj=None
	if os.path.getsize(filename)>0:
		if headerFlag==0:
		    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
		elif headerFlag==1:
		    fileDFObj = pd.read_table(filename, sep=sep, header=None, keep_default_na=False, na_values=[''])
		else:
		    print("Unrecognized Header Flag")
	return fileDFObj


def commonIds(idSet1, idSet2):
    ##this function finds ids commong in both Series.
    return idSet1[idSet1.isin(idSet2)]

def exclusive(idSet1, idSet2):
    ##this function finds ids commong in both Series.
    excSet1=idSet1[~idSet1.isin(idSet2)]
    excSet2=idSet2[~idSet2.isin(idSet1)]
    return pd.concat([excSet1,excSet2])

def fastaFilter(seqIds, fastFile):
    ##This function takes two parameters, one is the list of the ids and the other is the fasta file name.
    fastahandle = open(fastFile, "rU")
    #records = list(SeqIO.parse(fastahandle, "fasta"))
    #record_dict = SeqIO.index(fastFile, "fasta")
    record_dict = SeqIO.to_dict(SeqIO.parse(fastFile, "fasta"))
    fileFastaIds = pd.Series(list(record_dict.keys()))
    filteredIds=fileFastaIds[~fileFastaIds.isin(seqIds)]
    return filteredIds

def fastaFilterWrite(seqIds, fastFile, outFile):
    ##This function takes two parameters, one is the list of the ids and the other is the fasta file name.
    #fastahandle = open(fastFile, "rU")
    outHandle = open(outFile, "w")
    #records = list(SeqIO.parse(fastahandle, "fasta"))
    #record_dict = SeqIO.index(fastFile, "fasta")
    record_dict = SeqIO.to_dict(SeqIO.parse(fastFile, "fasta"))
    #print(list(record_dict.keys()))
    fileFastaIds = pd.Series(list(record_dict.keys()))
    filteredIds=fileFastaIds[fileFastaIds.isin(seqIds)]
    for i in range(0, filteredIds.size):
        #print(i)
        #print(filteredIds.iloc[i])
        #record_dict[filteredIds[i]]
        outHandle.write(">"+record_dict[filteredIds.iloc[i]].description+"\n")
        outHandle.write(str(record_dict[filteredIds.iloc[i]].seq)+"\n")
    outHandle.close()

def getORFs(outFile):
    ##This code runs transdecoder on the transcripts that either did not map to the geonome during PASA or failed the alignment quality test
    ##for both the aligner, gmap and blat.
    os.system(" "+ outFile)

def merge(outFile, pasa):
    ##This code merge the transcripts that either did not map to the geonome during PASA or failed the alignment quality test
    ##for both the aligner, gmap and blat and the pasa assembled transcripts. Then it launches the bowtie to map the reads to these transcripts.
    os.system("python merge_fasta_file.py "+ pasa+" "+ outFile+" "+  pasa+".trinity.fasta")

def runBowTie2(outFile, pasa):
    ##This code merge the transcripts that either did not map to the geonome during PASA or failed the alignment quality test
    ##for both the aligner, gmap and blat and the pasa assembled transcripts. Then it launches the bowtie to map the reads to these transcripts.
    os.system("sh ../shell/calculate "+ pasa+" "+ outFile+" "+  pasa+".trinity.fasta")

def main2(blatF, gmapF, blatV, gmapV, fastaFile):
    ##This main function reads blat and the gmap failed bed files, extracts trinity ids and calls commonIds and fastaFilter.
	blat1=readFile(blatF,'\t',1)
	gmap1=readFile(gmapF,'\t',1)
	blat2=readFile(blatV,'\t',1)
	gmap2=readFile(gmapV,'\t',1)
	if blat1 is not None:
		blat1Ids=blat1[3].str.extract("ID=([^,]+)")
		if gmap1 is not None:
			gmap1Ids=gmap1[3].str.extract("ID=([^,]+)")
			common=commonIds(blat1Ids, gmap1Ids)
			failedIds=pd.concat([blat1Ids, gmap1Ids])
		else:
			common=blat1Ids
			failedIds=blat1Ids
	else:
		if gmap1 is not None:
			gmap1Ids=gmap1[3].str.extract("ID=([^,]+)")
			common=gmap1Ids
			failedIds=gmap1Ids
		else:
			print("ERROR")
	
	if blat2 is not None:
		blat2Ids=blat2[3].str.extract("ID=([^,]+)")
		if gmap2 is not None:
			gmap2Ids=gmap2[3].str.extract("ID=([^,]+)")
			validIds=pd.concat([blat2Ids, gmap2Ids])
		else:
			validIds=blat2Ids
	else:
		if gmap2 is not None:
			gmap2Ids=gmap2[3].str.extract("ID=([^,]+)")
			validIds=gmap2Ids
		else:
			print("ERROR")

    #common=commonIds(blatIds, gmapIds)
	allIds=pd.concat([failedIds, validIds])
	unqAllIds=allIds.unique()
	notMapped=fastaFilter(unqAllIds, fastaFile)
	return pd.concat([common, notMapped])

def main(blatF, gmapF, trintyF, outF):
    ##This main function reads blat and the gmap failed bed files, extracts trinity ids and calls commonIds and fastaFilter.
    blat=readFile(blatF,'\t',1)
    gmap=readFile(gmapF,'\t',1)

    blatIds=blat[3].str.extract("ID=([^,]+)")
    gmapIds=gmap[3].str.extract("ID=([^,]+)")

    common=commonIds(blatIds, gmapIds)
    return common

parser = argparse.ArgumentParser(description='This python code find failed transcripts for both the gmap and blat aligner. These are possible virus transcripts. This code then extracts the trinity trascript and save in a file for further processing.')
parser.add_argument("-l", "--blatv", nargs=1, required=True, help="full path of the blat valid file", metavar="PATH")
parser.add_argument("-m", "--gmapv", nargs=1, required=True, help="full path of the gmap valid file", metavar="PATH")
parser.add_argument("-b", "--blatf", nargs=1, required=True, help="full path of the blat failed file", metavar="PATH")
parser.add_argument("-g", "--gmapf", nargs=1, required=True, help="full path of the gmap failed file", metavar="PATH")
parser.add_argument("-t", "--trinity", nargs=1, required=True, help="full path of to the trinity fasta", metavar="PATH")
parser.add_argument("-p", "--pasa", nargs=1, required=True, help="full path of to the PASA fasta", metavar="PATH")
parser.add_argument("-o", "--out", nargs=1, required=True, help="full path of to the output fasta file", metavar="PATH")

args = parser.parse_args()
print(args)

#blatF="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.failed_blat_alignments.bed"
#gmapF="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.failed_gmap_alignments.bed"
#blatV="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.valid_blat_alignments.bed"
#gmapV="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.valid_gmap_alignments.bed"
#trinityF="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adenovirus_trinity_assembled_transcripts.fasta"
#outFile="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/virus_and_novel.fasta"
virus_and_novel=main2(args.blatf[0], args.gmapf[0], args.blatv[0], args.gmapv[0], args.trinity[0])
fastaFilterWrite(virus_and_novel, args.trinity[0], args.out[0])
merge(args.out[0], args.pasa[0])


##adenoF="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/trinity/adenovirus_identified_orf_trinity_ids.csv"
##adeno=readFile(adenoF,"\t",1)
##idsList=adeno[0].str.split(",").tolist()
##ids=[item for sublist in idsList for item in sublist]

#python virusTranscriptExtractionAndORFPrediction.py -l /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.valid_blat_alignments.bed -m /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.valid_gmap_alignments.bed -b /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.failed_blat_alignments.bed -g /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.failed_gmap_alignments.bed -t /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adenovirus_trinity_assembled_transcripts.fasta -p /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/human_adeno_mydb_pasa.assemblies.fasta -o /data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/virus_and_novel2.fasta

#Namespace(blatf=['/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G69/G69.failed_blat_alignments.bed'], blatv=['/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G69/G69.valid_blat_alignments.bed'], gmapf=['/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G69/G69.failed_gmap_alignments.bed'], gmapv=['/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G69/G69.valid_gmap_alignments.bed'], out=['/data/SBCS-BessantLab/shyama/Data/Oliver/RSEM/G69/G69_missing.fasta'], pasa=['/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G69/G69.assemblies.fasta'], trinity=['/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/G69/G69.Trinity.fasta'])



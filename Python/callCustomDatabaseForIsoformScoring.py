##This code calls customDatabaseForIsoformScoring.py for each datasets.

import glob
import re
import os
##this code calles sample statistics for a given dataset.
##Mosquito
print("Mosquito")
orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PASA/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/identification/PASA/"
PITDB="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/PITDB/"
MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/MS/SliceAll.mgf"

s="aedes"
#command="python sampleStatistics.py --ORFs "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+"/"+s+".assemblies.fasta --proteins "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping+prt_filtered.csv --peptides "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping_filtered.csv --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsIdent "+PITDB+"transcripts/"+s+".assemblies.fasta.identified.fasta --ms "+MS+" --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"MosquitoSampleStat.tsv"
##ORFs Count	Transcript	Peptide	ORFs Count	Transcript	Peptide	ORFs Count	Transcript	Peptide	ORFs Count	Transcript	Peptide	ORFs Count	Transcript	Peptide	ORFs Count	Transcript	Peptide
command2="python sampleVariationStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"MosquitoSampleVarStat.tsv"
command3="python sampleIsoformStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf --vcfIso "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf >> "+PITDB+"MosquitoSampleIsoStat.tsv"

#os.system(command)
os.system(command2)
os.system(command3)

##Bat Nelson Bay
print("bat")
trinityDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/Trinity/"
orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PASA/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/"
PITDB="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/PITDB/"
MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/P.mgf"
for f in glob.glob(trinityDir+"*.Trinity.fasta"):
    s=re.sub("\.Trinity\.fasta","",os.path.basename(f))
    #command="python sampleStatistics.py --ORFs "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+"/"+s+".assemblies.fasta --proteins "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping+prt_filtered.csv --peptides "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping_filtered.csv --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsIdent "+PITDB+"transcripts/"+s+".assemblies.fasta.identified.fasta --ms "+MS+" --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"BatNelsonBaySampleStat.tsv"
    command2="python sampleVariationStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"BatNelsonBaySampleVarStat.tsv"
    command3="python sampleIsoformStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf --vcfIso "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf >> "+PITDB+"BatNelsonBaySampleIsoStat.tsv"
    if os.path.isfile(identified+s+"/"+s+"+fdr+th+grouping+prt.csv"):
        os.system(command2)
        os.system(command3)
        #os.system(command)


##Mouse nelson bay
print("Mouse")
trinityDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/Trinity/"
orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PASA/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/identification/"
PITDB="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/PITDB/"
MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/L.mgf"
for f in glob.glob(trinityDir+"*.Trinity.fasta"):
    s=re.sub("\.Trinity\.fasta","",os.path.basename(f))
	#print(s)
	#command="python sampleStatistics.py --ORFs "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+"/"+s+".assemblies.fasta --proteins "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping+prt_filtered.csv --peptides "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping_filtered.csv --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsIdent "+PITDB+"transcripts/"+s+".assemblies.fasta.identified.fasta --ms "+MS+" --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"MouseNelsonBaySampleStat.tsv"
    command2="python sampleVariationStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"MouseNelsonBaySampleVarStat.tsv"
    command3="python sampleIsoformStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf --vcfIso "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf >> "+PITDB+"MouseNelsonBaySampleIsoStat.tsv"
    if os.path.isfile(identified+s+"/"+s+"+fdr+th+grouping+prt.csv"):
		#print(command)
		#print("\n")
        #os.system(command)
        os.system(command2)
        os.system(command3)

##Human adeno-virus

'''
orfDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/"
transDir="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PASA/human_adeno_data/"
identified="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA/"
PITDB="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/PITDB/"
MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/MS/DM_from_raw.mgf"

s1="human_adeno_mydb_pasa"
s2="human_adeno"
#command="python sampleStatistics.py --ORFs "+orfDir+s1+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s1+".assemblies.fasta --proteins "+PITDB+"PSMs-Peptides-ORFs/"+s2+"+fdr+th+grouping+prt_filtered.csv --peptides "+PITDB+"PSMs-Peptides-ORFs/"+s2+"+fdr+th+grouping_filtered.csv --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s2+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsIdent "+PITDB+"transcripts/"+s2+".assemblies.fasta.identified.fasta --ms "+MS+" --annot "+PITDB+"Annotation/"+s1+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s2+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"HumanAdenoSampleStat.tsv"
command2="python sampleVariationStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s2+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s1+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s2+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"HumanAdenoSampleVarStat.tsv"
command2="python sampleVariationStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s2+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s1+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s2+".assemblies.fasta.transdecoder.pep_pepEvd.vcf --vcfIso "+PITDB+"Annotation/"+s1+".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf >> "+PITDB+"HumanAdenoSampleIsoStat.tsv"

#print(command)
#os.system(command)
#os.system(command2)
os.system(command3)
'''

##Oliver's data
print("Oliver")
msSamp={"G10":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S13.RAW.-1.mgf",
		"G11":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S1.RAW.-1.mgf",
		"G15":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S2.RAW.-1.mgf",
		"G17":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S3.RAW.-1.mgf",
		"G29a":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S4.RAW.-1.mgf",
		"G30":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S5.RAW.-1.mgf",
		"G33a":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S6.RAW.-1.mgf",
		"G36":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S14.RAW.-1.mgf",
		"G42":"3eab50639f039213___b1368p92_ECM_S26.RAW.-1.mgf",
		"G43":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S7.RAW.-1.mgf",
		"G54":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S8.RAW.-1.mgf",
		"G57":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S15.RAW.-1.mgf",
		"G58":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S9.RAW.-1.mgf",
		"G67":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S16.RAW.-1.mgf",
		"G69":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S17.RAW.-1.mgf",
		"G72":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S18.RAW.-1.mgf",
		"G75":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S19.RAW.-1.mgf",
		"G76":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S20.RAW.-1.mgf",
		"G79":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S22.RAW.-1.mgf",
		"G81":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S23.RAW.-1.mgf",
		"G82":"E_Vinni_Vinni_projects_B1368p20_b1368p20_sample82.mgf",
		"G85":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S10.RAW.-1.mgf",
		"G88":"3eab50639f039213___b1368p92_ECM_S29.RAW.-1.mgf",
		"G92":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S11.RAW.-1.mgf",
		"G93":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S12.RAW.-1.mgf",
		"G94":"3eab50639f039213___b1368p92_ECM_S27.RAW.-1.mgf",
		"G102":"E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S24.RAW.-1.mgf",
		"G124":"3eab50639f039213___b1368p92_ECM_S32.RAW.-1.mgf",
		"G138":"3eab50639f039213___b1368p92_ECM_S30.RAW.-1.mgf" }

MS="/data/SBCS-BessantLab/shyama/Data/Oliver/MS/"

trinityDir="/data/SBCS-BessantLab/shyama/Data/Oliver/Trinity/"
orfDir="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/"
transDir="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/"
identified="/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/"
PITDB="/data/SBCS-BessantLab/shyama/Data/Oliver/PITDB/"
#processedSamples=['G10','G11','G15','G17','G29a','G75','G102','G138']

for f in glob.glob(trinityDir+"*.Trinity.fasta"):
    print(f)
    s=re.sub("\.Trinity\.fasta","",os.path.basename(f))
    if os.path.exists(identified+s+"/"+s+"+fdr+th+grouping+prt.csv"):
        #print("\n"+s)
        #command="python sampleStatistics.py --ORFs "+orfDir+s+"/"+s+".assemblies.fasta.transdecoder.pep --transcripts "+transDir+s+"/"+s+".assemblies.fasta --proteins "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping+prt_filtered.csv --peptides "+PITDB+"PSMs-Peptides-ORFs/"+s+"+fdr+th+grouping_filtered.csv --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --transcriptsIdent "+PITDB+"transcripts/"+s+".assemblies.fasta.identified.fasta --ms "+MS+msSamp[s]+" --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"OliverSampleStat.tsv"
        command2="python sampleVariationStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf >> "+PITDB+"OliverSampleVarStat.tsv"
        command3="python sampleIsoformStats.py --ORFsIdent "+PITDB+"AminoAcids-or-ORFs-orTGEs/"+s+".assemblies.fasta.transdecoder.pep.identified.fasta --annot "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv --vcf "+PITDB+"Variations-proVCF/"+s+".assemblies.fasta.transdecoder.pep_pepEvd.vcf --vcfIso "+PITDB+"Annotation/"+s+".assemblies.fasta.transdecoder.pep_details_annotation.csv_isoform_pep.vcf>> "+PITDB+"OliverSampleIsoStat.tsv"
        
		#os.system(command)
        os.system(command2)
        os.system(command3)

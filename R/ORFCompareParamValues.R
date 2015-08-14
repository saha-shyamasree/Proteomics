### These are the parameter values of ORFCompare.R

#### Common values
upper=0.000000000000000000000000000001
pepThreshold=1
rev=1
peptide=1

####################################################  TrinityV5 Vs PASA ################################################
f2="pasa_assemblyV1+fdr+th+grouping+prt.csv"
f1="trinity_PITORF+fdr+th+grouping+prtV5.csv"
d2="D:/data/Results/Human-Adeno/Identification/PASA/sORF/"
d1="D:/data/Results/Human-Adeno/Identification/sORF/"
f4="human_adeno_mydb_pasa.assemblies_ORFsV1.csv"
f3="trinityV5Match.csv"
d4="D:/data/blast/blastCSV/PASA/Human-Adeno/"
d3="D:/data/blast/blastCSV/sORF/"
Mat1Name="TrinityV5"
Mat2Name="PASA"
outDir="D:/data/PASA/human_adeno/Stats/"
p2=lengthComparison(f1,f2,f3,f4,d1,d2,d3,d4,rev,peptide,pepThreshold,upper, Mat1Name, Mat2Name, outDir)


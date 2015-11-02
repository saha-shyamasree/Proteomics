## This Python code create the confidence Matrix.

prots="D:/data/Results/Human-Adeno/Identification/PASA/sORF/pasa_assemblyV1+fdr+th+grouping+prt.csv"
PSMs="D:/data/Results/Human-Adeno/Identification/PASA/sORF/pasa_assemblyV1+fdr+th+grouping.csv"
knownPro="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsV7.csv"
knownProSAP="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_knownProteinsSAPsV7.csv"
iso="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsV7.csv"
isoSAP="D:/data/blast/blastCSV/PASA/Human-Adeno/human_adeno_mydb_pasa.assemblies_ORFs_IsoformsSAPsV7.csv"
nonProt="D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonCodingPrt.tsv"
#nonPSM="D:/data/Results/Human-Adeno/Identification/PASA/sORF/bioTypeClassification/nonCodingPSM.tsv"
transcriptsFile=""
orfsFile=""
def readFile(filename, sep):
    fileDFObj = pd.read_table(filename, sep=sep, keep_default_na=False, na_values=[''])
    return fileDFObj

def pepThresholding(prots,pepTh):
    print("1. prots dim:"+str(prots.shape))
    prots=prots[~prots['protein accession'].str.contains("_REVERSED")]
    print("2. prots dim:"+str(prots.shape))
    prots=prots[prots['distinct peptide sequences']>pepTh]
    print("3. prots dim:"+str(prots.shape))
    return prots

def confidenceMatrix(confidenceFile, proteins, knownProt, knownSapProt, iso, isoSap, nonProt, variations, pepVariations, transcripts, orfs):
    for i in range(proteins):
        ## what to use? protein accession or description
        prot=proteins.iloc[i]['description']
        transcriptSeq=getTranscriptSeq(prot,transcripts)
        orfSeq=getORFSeq(prot,orfs)
        varCount=0
        pepVarCount=0
        protSimilarityClass=""
        biotype=""
        if checkProteinClass(knownProt,prot):
            protSimilarityClass="Known"
        elif checkProteinClass(knownSapProt,prot):
            protSimilarityClass="KnownVar"
            varCount=countVariations(prot, variations)
            pepVarCount=countVariations(prot, pepVariations)
        elif checkProteinClass(iso,prot):
            protSimilarityClass="Isoform"
        elif checkProteinClass(isoSap,prot):
            protSimilarityClass="IsoformVar"
            varCount=countVariations(prot, variations)
            pepVarCount=countVariations(prot, pepVariations)
        else:
            protSimilarityClass="Novel"
        
        writeToConfidenceFile(confidenceFile,proteins.iloc[i],orfSeq, transcriptSeq, protSimilarityClass,biotype, varCount, pepVarCount)
        
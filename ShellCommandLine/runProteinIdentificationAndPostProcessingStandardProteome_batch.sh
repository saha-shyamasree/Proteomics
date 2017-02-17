##This shell script call runProteinIdentificationAnsPostProcessing_cluster.py code for a project with multiple/single samples.
##This code has to be run from runProteinIdentificationAndPostProcessing_cluster.py directory
code='/data/home/btw796/Code2/Proteomics/Python'
<<COMMENT
## Oliver - human overian cancer (ECM)
standardPath="/data/SBCS-BessantLab/shyama/Data/Oliver/fasta/"
standardFile="uniprot-proteomeUP000005640.fasta"
ARRAY=( "G10:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S13.RAW.-1.mgf"
    "G11:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S1.RAW.-1.mgf"
    "G15:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S2.RAW.-1.mgf"
    "G17:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S3.RAW.-1.mgf"
    "G29a:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S4.RAW.-1.mgf"
    "G30:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S5.RAW.-1.mgf"
	"G33a:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S6.RAW.-1.mgf"
	"G36:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S14.RAW.-1.mgf"
	"G42:3eab50639f039213___b1368p92_ECM_S26.RAW.-1.mgf"
	"G43:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S7.RAW.-1.mgf"
	"G54:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S8.RAW.-1.mgf"
	"G57:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S15.RAW.-1.mgf"
	"G58:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S9.RAW.-1.mgf"
	"G67:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S16.RAW.-1.mgf"
	"G69:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S17.RAW.-1.mgf"
	"G72:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S18.RAW.-1.mgf"
    "G75:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S19.RAW.-1.mgf"
	"G76:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S20.RAW.-1.mgf"
	"G79:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S22.RAW.-1.mgf"
	"G81:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S23.RAW.-1.mgf"
	"G82:E_Vinni_Vinni_projects_B1368p20_b1368p20_sample82.mgf"
	"G85:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S10.RAW.-1.mgf"
	"G88:3eab50639f039213___b1368p92_ECM_S29.RAW.-1.mgf"
	"G92:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S11.RAW.-1.mgf"
	"G93:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S12.RAW.-1.mgf"
	"G94:3eab50639f039213___b1368p92_ECM_S27.RAW.-1.mgf"
    "G102:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S24.RAW.-1.mgf"
	"G124:3eab50639f039213___b1368p92_ECM_S32.RAW.-1.mgf"
	"G138:3eab50639f039213___b1368p92_ECM_S30.RAW.-1.mgf" )

MS="/data/SBCS-BessantLab/shyama/Data/Oliver/MS/"
ident="/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/"
modification="/data/SBCS-BessantLab/shyama/Data/Oliver/MS/modifications.txt"
COMMENT

## Bristol
##Mouse - Nelson Bay
#standardPath="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/fasta/"
#standardFile="MouseNelsonBayUniprot.fasta"
#ARRAY=( "mouseNelsonBay:L.mgf")
#MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/MS/"
#ident="/data/SBCS-BessantLab/shyama/Data/Bristol/Mouse/identification/"
#modification="/data/SBCS-BessantLab/shyama/Data/Bristol/modifications.txt"

##Mosquito
standardPath="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/fasta/"
standardFile="Aedes_aegypti_uniprot-proteomeUP000008820.fasta"
ARRAY=( "mosquito:SliceAll.mgf")
MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/MS/"
ident="/data/SBCS-BessantLab/shyama/Data/Bristol/Mosquito/identification/"
modification="/data/SBCS-BessantLab/shyama/Data/Bristol/modifications_mosquito.txt"

##Bat - Nelson Bay
#standardPath="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/fasta/"
#standardFile="Palecto_nelsonBay.fasta"
#ARRAY=( "batNelson:P.mgf")
#MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/MS/"
#ident="/data/SBCS-BessantLab/shyama/Data/Bristol/Bat/NelsonBay/identification/"
#modification="/data/SBCS-BessantLab/shyama/Data/Bristol/modifications.txt"
## Human - adeno
#standard="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/fasta/HumanAdenoPapilloma.fasta"
#ARRAY=( "human_adeno_mydb_pasa:DM_from_raw.mgf")
#
#MS="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/MS/"
#ident="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/identification/PASA"
#modification="/data/SBCS-BessantLab/shyama/Data/Bristol/Human/adenovirus/MS/modifications.txt"
for arr in "${ARRAY[@]}" ; do
    sample=${arr%%:*}
    mFile=${arr#*:}
    if [ ! -d $ident$sample ]; then
        ##comment following statement when an experiment has one sample
    	mkdir $ident$sample
        #echo "not required"
    fi
    #echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i $MS$mFile -o $ident$sample/$sample.mzid -d $PASA$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Oliver/MS/modifications.txt > $ident$sample/run_$sample.cmds.sh"
    ##Oliver
    ##Check the standard file , cause all of them are using same file, that may create problem. Because MSGF+ index/or somthing the fasta file.
    cp $standardPath$standardFile $ident$sample/$standardFile
    python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i $MS$mFile -o $ident$sample/$sample.standard.mzid -d $ident$sample/$standardFile -m $modification > $ident$sample/run_$sample.standard.cmds.sh
    ##Run this when an experiment has one sample, e.g. human adeno, mosquito
    ##Human adeno, modification files are different
    #python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i $MS$mFile -o $ident/$sample.standard.mzid -d $standard -m $modification > $ident/run_$sample.standard.cmds.sh
    
    printf "%s sample %s\n" "$sample" "$mFile"
    qsub $ident$sample/run_$sample.standard.cmds.sh
    ##Run this when an experiment has one sample, e.g. human adeno, mosquito
    #qsub $ident/run_$sample.standard.cmds.sh
    #break
done


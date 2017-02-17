##This shell script call runProteinIdentificationAnsPostProcessing_cluster.py code for a project with multiple samples.
##This code has to be run from runProteinIdentificationAndPostProcessing_cluster.py directory
code='/data/home/btw796/Code2/Proteomics/Python'
## Oliver - human overian cancer (ECM)
ARRAY=( "G30:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S5.RAW.-1.mgf"
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
	"G76:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S20.RAW.-1.mgf"
	"G79:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S22.RAW.-1.mgf"
	"G81:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S23.RAW.-1.mgf"
	"G82:E_Vinni_Vinni_projects_B1368p20_b1368p20_sample82.mgf"
	"G85:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S10.RAW.-1.mgf"
	"G88:3eab50639f039213___b1368p92_ECM_S29.RAW.-1.mgf"
	"G92:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S11.RAW.-1.mgf"
	"G93:E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S12.RAW.-1.mgf"
	"G94:3eab50639f039213___b1368p92_ECM_S27.RAW.-1.mgf"
	"G124:3eab50639f039213___b1368p92_ECM_S32.RAW.-1.mgf"
	"G138:3eab50639f039213___b1368p92_ECM_S30.RAW.-1.mgf" )

MS="/data/SBCS-BessantLab/shyama/Data/Oliver/MS/"
ident="/data/SBCS-BessantLab/shyama/Data/Oliver/Identification/"
PASA="/data/SBCS-BessantLab/shyama/Data/Oliver/PASA/"

for arr in "${ARRAY[@]}" ; do
    sample=${arr%%:*}
    mFile=${arr#*:}
    if [ ! -d $ident$sample ]; then
    	mkdir $ident$sample
    fi
    #echo "python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i $MS$mFile -o $ident$sample/$sample.mzid -d $PASA$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Oliver/MS/modifications.txt > $ident$sample/run_$sample.cmds.sh"
    #python3.4 runProteinIdentificationAndPostProcessing_cluster.py -i $MS$mFile -o $ident$sample/$sample.mzid -d $PASA$sample/$sample.assemblies.fasta.transdecoder.pep -m /data/SBCS-BessantLab/shyama/Data/Oliver/MS/modifications.txt > $ident$sample/run_$sample.cmds.sh
    
    printf "%s sample %s\n" "$sample" "$mFile"
    qsub $ident$sample/run_$sample.cmds.sh
    #break
done


################################## MSConvert  ##########################################

############################### MSGF+ ######################################
#Update 09/06/2015
# I should allow at least 2 maps for a spectrum so that I can compre the delta score.

##############Uniprot Second run with adenovirus falimy and including Human Papiloma Virus #####################
java -Xmx12500M -jar MSGFPlus.jar -s J:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\uniprot_Human_Mastadenovirus_human_papiloma_virus_types_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2.mzid -mod J:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms10024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms5024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms5024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr+th.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
#java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv C:\Data\HUMAN\mgf\DM_from_raw_uni_trinity+fdr+th+grouping.mzid C:\Data\HUMAN\mgf\DM_from_raw_uni_trinity+fdr+th+grouping_prt_grp.csv -exportType  exportProteinGroups  -compress false
java -Xms5024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms5024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot2+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

##uniprot and trinity
java -Xmx8500M -jar MSGFPlus.jar -s C:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\Trinity\trinity_uniprotV2_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms1024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr+th.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
#java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv C:\Data\HUMAN\mgf\DM_from_raw_uni_trinity+fdr+th+grouping.mzid C:\Data\HUMAN\mgf\DM_from_raw_uni_trinity+fdr+th+grouping_prt_grp.csv -exportType  exportProteinGroups  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false
##uniprot and trinity without duplicates

java -Xmx8500M -jar MSGFPlus.jar -s C:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\Trinity\trinity_uniprot_no_duplicateV2_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms1024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr+th.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity_uniprot_no_duplicate+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false
##trinity

java -Xmx8500M -jar MSGFPlus.jar -s C:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\Trinity\trinityV2_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\DM_from_raw_trinity.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms1024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\DM_from_raw_trinity.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr+th.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_trinity+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

### PITORFAll
Supplementary_dataset_1_-_human_trinity_assembled_transcripts-ORF_concatenated_target_decoy

java -Xmx8500M -jar MSGFPlus.jar -s C:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\Trinity\Supplementary_dataset_1_-_human_trinity_assembled_transcripts-ORF_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\trinity_PITORF.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms1024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\trinity_PITORF.mzid G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr.mzid G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr+th.mzid G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\trinity_PITORF+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false
#on laptop
java -Xms1024m -jar mzidentml-lib.jar Threshold D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr.mzid D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr+th2.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold false
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr+th2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr+th+grouping2.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr+th+grouping2.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\trinity\trinity_PITORF+fdr+th+grouping+prt2.csv -exportType  exportProteinsOnly  -compress false


###SORF shortest mainORF - 30bp, shortest uORF 24bp, longest 5` uORF - 200bp
java -Xmx12500M -jar MSGFPlus.jar -s J:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\Trinity\trinityV3_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF.mzid -mod J:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms8024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms8024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms8024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms8024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms8024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

##sORFs, shortest mainORF - 30bp, shortest uORF 24bp, longest 5` uORF - 300bp
java -Xmx12500M -jar MSGFPlus.jar -s J:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\Trinity\trinityV4_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORFV4.mzid -mod J:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms8024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORFV4.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdrV4.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms8024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdrV4.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+thV4.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms8024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+thV4.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV4.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms8024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV4.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV4.csv -exportType  exportPSMs  -compress false
java -Xms8024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV4.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+grouping+prtV4.csv -exportType  exportProteinsOnly  -compress false


##sORFs, shortest mainORF - 30bp, shortest uORF 24bp, longest 5` uORF - 600bp
java -Xmx12500M -jar MSGFPlus.jar -s J:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\Trinity\trinityV5_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORFV5.mzid -mod J:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms8024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORFV5.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdrV5.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms8024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdrV5.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+thV5.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms8024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+thV5.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV5.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms8024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV5.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV5.csv -exportType  exportPSMs  -compress false
java -Xms8024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+groupingV5.mzid G:\SearchEngine_output\MSGF+\HumanAdeno\sORF\trinity_PITORF+fdr+th+grouping+prtV5.csv -exportType  exportProteinsOnly  -compress false

##uniprot

java -Xmx8500M -jar MSGFPlus.jar -s C:\Data\HUMAN\mgf\DM_from_raw.mgf -d K:\E\Data\HUMAN\database\human_adenovirus_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms1024m -jar mzidentml-lib.jar FalseDiscoveryRate C:\Data\HUMAN\mgf\DM_from_raw_uniprot.mzid C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr.mzid C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true

java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr+th.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot+fdr+th+prt.csv -exportType  exportProteinsOnly  -compress false

java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr+th.mzid C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr+th+grouping.mzid C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv C:\Data\HUMAN\mgf\DM_from_raw_uniprot+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\DM_from_raw_uniprot+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false


#on laptop, not removing PSMs below threshold
java -Xms1024m -jar mzidentml-lib.jar Threshold D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr+th2.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold false
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr+th2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr+th+grouping2.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr+th+grouping2.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-trinity\DM_from_raw_uniprot+fdr+th+grouping+prt2.csv -exportType  exportProteinsOnly  -compress false


##Cufflinks (genome guided, main-sORFs included)
java -Xmx10500M -jar MSGFPlus.jar -s J:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\cufflinks\transcripts-ORF_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\cufflinks_PITORF.mzid -mod J:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms5024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\cufflinks_PITORF.mzid G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr.mzid G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr+th.mzid G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\cufflinks_PITORF+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

##Cufflinks (genome guided, main ORFs)
java -Xmx10500M -jar MSGFPlus.jar -s J:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\cufflinks\transcripts_main-ORF_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF.mzid -mod J:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms5024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF.mzid G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr.mzid G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr+th.mzid G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\cufflinks_main_PITORF+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

#on laptop, not removing PSMs below threshold
java -Xms1024m -jar mzidentml-lib.jar Threshold D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr.mzid D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr+th2.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold false
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr+th2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr+th+grouping2.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr+th+grouping2.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\cufflinks\cufflinks_main_PITORF+fdr+th+grouping+prt2.csv -exportType  exportProteinsOnly  -compress false

##Uniprot, genome guided
### Only human uniprot ###

java -Xmx12500M -jar MSGFPlus.jar -s J:\Data\HUMAN\mgf\DM_from_raw.mgf -d E:\Data\HUMAN\database\HUMAN_concatenated_target_decoy.fasta -o G:\SearchEngine_output\MSGF+\human_uniprot.mzid -mod J:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8 
java -Xms10024m -jar mzidentml-lib.jar FalseDiscoveryRate G:\SearchEngine_output\MSGF+\human_uniprot.mzid G:\SearchEngine_output\MSGF+\human_uniprot+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms5024m -jar mzidentml-lib.jar Threshold G:\SearchEngine_output\MSGF+\human_uniprot+fdr.mzid G:\SearchEngine_output\MSGF+\human_uniprot+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms5024m -jar mzidentml-lib.jar ProteoGrouper G:\SearchEngine_output\MSGF+\human_uniprot+fdr+th.mzid G:\SearchEngine_output\MSGF+\human_uniprot+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
#java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv C:\Data\HUMAN\mgf\DM_from_raw_uni_trinity+fdr+th+grouping.mzid C:\Data\HUMAN\mgf\DM_from_raw_uni_trinity+fdr+th+grouping_prt_grp.csv -exportType  exportProteinGroups  -compress false
java -Xms5024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\human_uniprot+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\human_uniprot+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms5024m -jar mzidentml-lib.jar Mzid2Csv G:\SearchEngine_output\MSGF+\human_uniprot+fdr+th+grouping.mzid G:\SearchEngine_output\MSGF+\human_uniprot+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

#on Laptop, not removing PSMs below threshold

java -Xms1024m -jar mzidentml-lib.jar Threshold D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr+th2.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold false
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr+th2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr+th+grouping2.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr+th+grouping2.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr+th+grouping2.mzid D:\data\Results\Human-Adeno\GIOPaperResults\Uniprot-Cufflinks\human_uniprot+fdr+th+grouping+prt2.csv -exportType  exportProteinsOnly  -compress false

############################################## X!Tandem ##########################################

##Ran X!Tandem from SearchGUI on my local machines, parameter files are saved
##convert Tandem output to mzid format so that I can analysis the data using same pipeline and also merge the results with other search engines results.
java -Xms1024m -jar mzidentml-lib.jar Tandem2mzid G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_uniprot.t.xml G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_uniprot.mzid -decoyRegex _REVERSED -databaseFileFormatID "MS:1001348" -massSpecFileFormatID "MS:1001062" -idsStartAtZero false -compress false
java -Xms1024m -jar mzidentml-lib.jar Tandem2mzid G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_trinity.t.xml G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_trinity.mzid -decoyRegex _REVERSED -databaseFileFormatID "MS:1001348" -massSpecFileFormatID "MS:1001062" -idsStartAtZero false -compress false
java -Xms1024m -jar mzidentml-lib.jar Tandem2mzid G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_uni_trinity.t.xml G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_trinity_uniprot.mzid -decoyRegex _REVERSED -databaseFileFormatID "MS:1001348" -massSpecFileFormatID "MS:1001062" -idsStartAtZero false -compress false
java -Xms1024m -jar mzidentml-lib.jar Tandem2mzid G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_uni_trinity_no_duplicate.t.xml G:\SearchEngine_output\X!Tandem\PIT_data_human\DM_from_raw_trinity_uniprot_no_duplicate.mzid -decoyRegex _REVERSED -databaseFileFormatID "MS:1001348" -massSpecFileFormatID "MS:1001062" -idsStartAtZero false -compress false


##################################################### Bat data ################################################

I:\Bat\Data\Bat_HendraSILAC\Repeat
######################    MSConvert   ####################
msconvert -f I:\Bat\Data\Bat_HendraSILAC\Repeat\filenames.txt -o I:\Bat\Data\mgf --mgf --merge
################  MSGF+  #################
####  Trinity  ####
java -Xmx12000M -jar MSGFPlus.jar -s I:\Bat\Data\mgf\Slice.mgf -d I:\Bat\Data\Trinity_P_alecto-ORF_concatenated_target_decoy.fasta -o I:\Bat\SearchEngine\MSGF+\trinity_bat.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8
java -Xms3024m -jar mzidentml-lib.jar FalseDiscoveryRate I:\Bat\SearchEngine\MSGF+\trinity_bat.mzid I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms1024m -jar mzidentml-lib.jar Threshold I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr.mzid I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th.mzid I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th+grouping.mzid I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th+grouping.mzid I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th+grouping.mzid I:\Bat\SearchEngine\MSGF+\trinity_bat+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false
######### uniprot  ####### 
java -Xmx12000M -jar MSGFPlus.jar -s I:\Bat\Data\mgf\Slice.mgf -d I:\Bat\Data\uniprot_concatenated_target_decoy.fasta -o I:\Bat\SearchEngine\MSGF+\uniprot_bat.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8
java -Xms10000m -jar mzidentml-lib.jar FalseDiscoveryRate I:\Bat\SearchEngine\MSGF+\uniprot_bat.mzid I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms3024m -jar mzidentml-lib.jar Threshold I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr.mzid I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms1024m -jar mzidentml-lib.jar ProteoGrouper I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr+th.mzid I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr+th+grouping.mzid I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms1024m -jar mzidentml-lib.jar Mzid2Csv I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr+th+grouping.mzid I:\Bat\SearchEngine\MSGF+\uniprot_bat+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false



#########################################################       Human with Hendra         ##########################################################
msconvert -f I:\Human-Hendra\Data\filenames.txt -o I:\Human-Hendra\Data\mgf --mgf --merge

######  MSGF+  #####
#####trinity#####
java -Xmx12000M -jar MSGFPlus.jar -s I:\Human-Hendra\Data\mgf\Slice.mgf -d I:\Human-Hendra\Trinity_HeV293-ORF_concatenated_target_decoy.fasta -o I:\Human-Hendra\SearchEngine\MSGF+\trinity.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8
java -Xms10024m -jar mzidentml-lib.jar FalseDiscoveryRate I:\Human-Hendra\SearchEngine\MSGF+\trinity.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms8024m -jar mzidentml-lib.jar Threshold I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms4024m -jar mzidentml-lib.jar ProteoGrouper I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

#####Uniprot####
java -Xmx11000M -jar MSGFPlus.jar -s I:\Human-Hendra\Data\mgf\Slice.mgf -d I:\Human-Hendra\uniprot_concatenated_target_decoy.fasta -o I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8
java -Xms11024m -jar mzidentml-lib.jar FalseDiscoveryRate I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human.mzid I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms4024m -jar mzidentml-lib.jar Threshold I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr.mzid I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms4024m -jar mzidentml-lib.jar ProteoGrouper I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr+th.mzid I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\uniprot_human+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

############################################# Human Adeno, PASA assembly ORFs   ################################################
I:\Human-Adenovirus\ProteinIdent
java -Xmx12000M -jar MSGFPlus.jar -s /data/scratch/btw796/Data/HumanAdeno/MS/DM_from_raw.mgf -d  /data/scratch/btw796/Data/HumanAdeno/fasta/human_adeno_mydb_pasa.assemblies_ORFs_concatenated_target_decoy.fasta -o I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1.mzid -mod /data/scratch/btw796/Data/HumanAdeno/MS/modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8
java -Xmx10024m -jar mzidentml-lib.jar FalseDiscoveryRate I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1.mzid I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xmx8024m -jar mzidentml-lib.jar Threshold I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr.mzid I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xmx4024m -jar mzidentml-lib.jar ProteoGrouper I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr+th.mzid I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xmx4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr+th+grouping.mzid I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xmx4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr+th+grouping.mzid I:\Human-Adenovirus\ProteinIdent\pasa_assemblyV1+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false
## This python code run MSGF+ and mzidenML-lib i.e. protein identification and post-process the result

import argparse
import subprocess as sp

parser = argparse.ArgumentParser(description='This python code run MSGF+ and mzidenML-lib i.e. protein identification and post-process the result')
parser.add_argument("-i", "--input", nargs=1, required=True, help="full path of mgf file", metavar="PATH")
parser.add_argument("-o", "--output", nargs=1, required=True, help="full path of of output file", metavar="PATH")
parser.add_argument("-d", "--database", nargs=1, required=True, help="full path to the database file", metavar="PATH")
parser.add_argument("-m","--modification", nargs=1, required=True, help="full path to the modification file", metavar="PATH")
parser.add_argument("-t", "--tolerance", default="10ppm", help="mass tolerance")
#parser.add_option("--tda", help="whether to create decoy database", metavar="Integer")
parser.add_argument("--m_msgf", default=1, help="", metavar="Integer")
parser.add_argument("--inst", default=1,help="")
parser.add_argument("--minLength", default=8, help="minimum peptide length", metavar="LENGTH")
parser.add_argument("--msgf_path", default="",help="msgf+ jar file location", metavar="PATH")
#parser.add_option("--min3", help="the minimum length of 3' uORFs, only works when -uorf3output is set", metavar="Integer")

args = parser.parse_args()
print(args)

'''

    name or flags - Either a name or a list of option strings, e.g. foo or -f, --foo.
    action - The basic type of action to be taken when this argument is encountered at the command line.
    nargs - The number of command-line arguments that should be consumed.
    const - A constant value required by some action and nargs selections.
    default - The value produced if the argument is absent from the command line.
    type - The type to which the command-line argument should be converted.
    choices - A container of the allowable values for the argument.
    required - Whether or not the command-line option may be omitted (optionals only).
    help - A brief description of what the argument does.
    metavar - A name for the argument in usage messages.
    dest - The name of the attribute to be added to the object returned by parse_args().

'''

#subprocess.call(["perl /data/home/btw796/Code/Proteomics/Perl/orfall.pl", "-h"])
sp.check_output(["cd",args.msgf_path], shell=True, stderr=subprocess.STDOUT)
sp.check_output(["java", "-Xmx12000M", "-jar", "MSGFPlus.jar","D:\Code\Proteomics\Perl\orfall.pl",args.input,args.output])


'''
java -Xmx12000M -jar MSGFPlus.jar -s I:\Human-Hendra\Data\mgf\Slice.mgf -d I:\Human-Hendra\Trinity_HeV293-ORF_concatenated_target_decoy.fasta -o I:\Human-Hendra\SearchEngine\MSGF+\trinity.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8
java -Xms10024m -jar mzidentml-lib.jar FalseDiscoveryRate I:\Human-Hendra\SearchEngine\MSGF+\trinity.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms8024m -jar mzidentml-lib.jar Threshold I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms4024m -jar mzidentml-lib.jar ProteoGrouper I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false

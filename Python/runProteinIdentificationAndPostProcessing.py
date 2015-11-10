## This python code run MSGF+ and mzidenML-lib i.e. protein identification and post-process the result

import argparse
from subprocess import Popen, PIPE
#import subprocess as sp
import os

parser = argparse.ArgumentParser(description='This python code run MSGF+ and mzidenML-lib i.e. protein identification and post-process the result')
parser.add_argument("-i", "--input", nargs=1, required=True, help="full path of mgf file", metavar="PATH")
parser.add_argument("-o", "--output", nargs=1, required=True, help="full path of of output file", metavar="PATH")
parser.add_argument("-d", "--database", nargs=1, required=True, help="full path to the database file", metavar="PATH")
parser.add_argument("-m","--modification", nargs=1, required=True, help="full path to the modification file", metavar="PATH")
parser.add_argument("-t", "--tolerance", default=["10ppm"], help="mass tolerance")
#parser.add_option("--tda", help="whether to create decoy database", metavar="Integer")
parser.add_argument("--m_msgf", default=['1'], help="", metavar="Integer")
parser.add_argument("--inst", default=['1'],help="Instrument")
parser.add_argument("--minLength", default=['8'], help="minimum peptide length", metavar="LENGTH")
parser.add_argument("--tda", default=['1'],help="Decoy database search tda 0/1(default)")

parser.add_argument("--msgf_path", default=["C:/soft/MSGFPlus.20140210/"],help="msgf+ jar file location", metavar="PATH")
parser.add_argument("--mzident_path", default=["C:/soft/MzidLib-1.6/"],help="mzidentml-lib jarfile location", metavar="PATH")
#parser.add_option("--min3", help="the minimum length of 3' uORFs, only works when -uorf3output is set", metavar="Integer")

args = parser.parse_args()
print(args)

#contaminants="/data/SBCS-BessantLab/shyama/Data/Contaminants/crap.fasta"
contaminants="J:\Oliver\PASA\\apocrita\crap.fasta"

com1=" ".join(["python","merge_fasta_file.py", args.database[0], contaminants,  args.database[0]+".cont.fasta"])
#p=sp.call(["python3.4","merge_fasta_file.py", args.database[0], contaminants,  args.database[0]+".cont.fasta"])
#p.wait()
#outMsg = p.communicate()[0]
#print(p.returncode)
'''
os.system("python merge_fasta_file.py "+args.database[0]+" "+contaminants+" "+args.database[0]+".cont.fasta")
print("ran merge code")
'''

os.chdir(args.msgf_path[0])
print("current dir:")
print(os.getcwd())
command=" ".join(["java", "-Xmx12000M", "-jar", "MSGFPlus.jar","-s",args.input[0],"-d",args.database[0]+".cont.fasta","-o",args.output[0],"-mod",args.modification[0],"-t",args.tolerance[0],"-m",args.m_msgf[0],"-tda",args.tda[0],"-inst",args.inst[0],"-minLength",args.minLength[0]])
print(command) #"-Xmx5000m",
#os.system("java -XX:+PrintFlagsFinal -version | findstr HeapSize")
#os.system("java -verbose")
os.system(command)
#msgf=Popen(["java","-Xmx6100M", "-Xms6100M","-jar", "MSGFPlus.jar","-s",args.input[0],"-d",args.database[0]+".cont.fasta","-o",args.output[0],"-mod",args.modification[0],"-t",args.tolerance[0],"-m",args.m_msgf[0],"-tda",args.tda[0],"-inst",args.inst[0],"-minLength",args.minLength[0]])
#msgf.wait()

os.chdir(args.mzident_path[0])
print(os.getcwd())
mzFDR=" ".join(["java","-Xms1024m","-jar","mzidentml-lib.jar","FalseDiscoveryRate",args.database[0]+".mzid",args.database[0]+"+fdr.mzid","-decoyRegex","XXX_","-decoyValue","1","-cvTerm","\"MS:1002053\"","-betterScoresAreLower","true"])
print(mzFDR)
os.system(mzFDR)

mzTh=" ".join(["java","-Xmx8024m","-jar","mzidentml-lib.jar","Threshold",args.database[0]+"+fdr.mzid",args.database[0]+"+fdr+th.mzid","-isPSMThreshold","true","-cvAccessionForScoreThreshold","MS:1002355","-threshValue","0.01","-betterScoresAreLower","true","-deleteUnderThreshold","true"])
os.system(mzTh)

mzGrp=" ".join(["java","-Xmx4024m","-jar","mzidentml-lib.jar","ProteoGrouper",args.database[0]+"+fdr+th.mzid",args.database[0]+"+fdr+th+grouping.mzid","-cvAccForSIIScore","MS:1002355","-logTransScore","true","-requireSIIsToPassThreshold","true","-verboseOutput","false"])
os.system(mzGrp)

mzPep=" ".join(["java","-Xmx4024m","-jar","mzidentml-lib.jar","Mzid2Csv",args.database[0]+"+fdr+th+grouping.mzid",args.database[0]+"+fdr+th+grouping.csv","-exportType","exportPSMs","-compress","false"])
os.system(mzPep)

mzPrt=" ".join(["java","-Xmx4024m","-jar","mzidentml-lib.jar","Mzid2Csv",args.database[0]+"+fdr+th+grouping.mzid",args.database[0]+"+fdr+th+grouping+prt.csv","-exportType","exportProteinsOnly","-compress","false"])
os.system(mzPrt)

'''
java -Xmx12000M -jar MSGFPlus.jar -s I:\Human-Hendra\Data\mgf\Slice.mgf -d I:\Human-Hendra\Trinity_HeV293-ORF_concatenated_target_decoy.fasta -o I:\Human-Hendra\SearchEngine\MSGF+\trinity.mzid -mod C:\Data\HUMAN\mzML\modifications.txt -t 10ppm -m 1 -tda 0 -inst 1 -minLength 8
java -Xms10024m -jar mzidentml-lib.jar FalseDiscoveryRate I:\Human-Hendra\SearchEngine\MSGF+\trinity.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr.mzid -decoyRegex _REVERSED -decoyValue 1 -cvTerm "MS:1002053" -betterScoresAreLower true
java -Xms8024m -jar mzidentml-lib.jar Threshold I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1002355 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true
java -Xms4024m -jar mzidentml-lib.jar ProteoGrouper I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid -cvAccForSIIScore MS:1002355 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.csv -exportType  exportPSMs  -compress false
java -Xms4024m -jar mzidentml-lib.jar Mzid2Csv I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping.mzid I:\Human-Hendra\SearchEngine\MSGF+\trinity+fdr+th+grouping+prt.csv -exportType  exportProteinsOnly  -compress false
'''

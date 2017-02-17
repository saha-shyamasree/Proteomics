## this code calles PreliminaryProteinAnnotationForPITDB.py


import os
import argparse
import glob

parser = argparse.ArgumentParser(description='This python code calls  PreliminaryProteinAnnotationForPITDB.py.')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="full path of to the Location folder of the identified ORFs( possibly within Blast folder", metavar="PATH")
#parser.add_argument("-u", "--uniprot", nargs=1, required=True, help="full path of uniprot file", metavar="PATH")
args = parser.parse_args()
#print(args)
#onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
folders=glob.glob(args.blast[0]+"*")
#onlyfiles=glob.glob(args.blast[0]+"/*.assemblies.fasta.transdecoder.pep.identified.loc.csv")
#print(args.blast[0])
for f in folders:
    #print(f)
    sample=os.path.basename(f)
    #print(sample)
    #print(args.uniprot[0])
    #fBase=str(os.path.join(args.blast[0],f)).rstrip(".identified.loc.csv")
    #cmd="python PreliminaryProteinAnnotationForPITDB.py -b "+f+"/"+sample+".assemblies.fasta.transdecoder.pep.identified.loc.csv"
    cmd="python PreliminaryProteinAnnotationForPITDBV2.py -b "+f+"/"+sample+".assemblies.fasta.transdecoder.pep.identified.loc.csv"
    # -k "+f+"/"+sample+".assemblies.fasta.transdecoder.pep_known.csv -s "+f+"/"+sample+".assemblies.fasta.transdecoder.pep_knownVar.csv -v "+f+"/"+sample+".assemblies.fasta.transdecoder.pep.vcf -i "+f+"/"+sample+".assemblies.fasta.transdecoder.pep_iso.csv -j "+f+"/"+sample+".assemblies.fasta.transdecoder.pep_isoVar.csv"
    #print(cmd)
    os.system(cmd)
    #break

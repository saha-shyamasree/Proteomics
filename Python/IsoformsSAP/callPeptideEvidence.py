## this code calls peptideEvidence.py for different samples

import os
import argparse
import glob

parser = argparse.ArgumentParser(description='This python code run IdentifyProteinIsoformSAP.py')
parser.add_argument("-p", "--psm", nargs=1, required=True, help="full path of blast folder", metavar="PATH")
parser.add_argument("-v", "--vcf", nargs=1, required=True, help="full path of blast folder", metavar="PATH")
#parser.add_argument("-u", "--uniprot", nargs=1, required=True, help="full path of uniprot file", metavar="PATH")
args = parser.parse_args()
print(args)
#onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
onlyfiles=glob.glob(args.psm[0]+"/*.assemblies.fasta.transdecoder.pep+fdr+th+grouping.csv")
print(args.psm[0])
for f in onlyfiles:
    print("F:"+f)
    #print(args.uniprot[0])
    fBase=f.rstrip(".csv")
    #print("fBase"+fBase)
    sample=os.path.basename(f).split(".",maxsplit=1)[0]
    print("sample:"+str(sample))
    cmd="python D:\Code\Proteomics\Python\IsoformsSAP\peptideEvidence.py -p "+f+" -s "+fBase+"_Variation.csv -v "+args.vcf[0]+"/"+sample+".assemblies.fasta.transdecoder.pep.vcf -o "+args.vcf[0]+"/"+sample+".assemblies.fasta.transdecoder.pep_pepEvd.vcf"
    print(cmd)
    #print("\n")
    os.system(cmd)
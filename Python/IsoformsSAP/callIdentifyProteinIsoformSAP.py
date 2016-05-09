## this code calls IdentifyProteinIsoformSAP.py for different samples

import os
import argparse
import glob

parser = argparse.ArgumentParser(description='This python code run IdentifyProteinIsoformSAP.py')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="full path of blast folder", metavar="PATH")
#parser.add_argument("-u", "--uniprot", nargs=1, required=True, help="full path of uniprot file", metavar="PATH")
args = parser.parse_args()
print(args)
#onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
onlyfiles=glob.glob(args.blast[0]+"/*.assemblies.fasta.transdecoder.pep.csv")
print(args.blast[0])
for f in onlyfiles:
    print(f)
    #print(args.uniprot[0])
    fBase=str(os.path.join(args.blast[0],f)).rstrip(".csv")
    cmd="python D:\Code\Proteomics\Python\IsoformsSAP\IdentifyProteinIsoformSAP.py -b "+str(os.path.join(args.blast[0],f))+" -k "+fBase+"_known.csv"+" -s "+fBase+"_knownVar.csv -v "+fBase+".vcf -i "+fBase+"_iso.csv -j "+fBase+"_isoVar.csv"
    print(cmd)
    os.system(cmd)
    #break

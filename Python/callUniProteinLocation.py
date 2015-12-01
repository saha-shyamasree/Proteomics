## this code calls UniProteinLocation.py for different samples

import os
import argparse

parser = argparse.ArgumentParser(description='This python code run UniProteinLocation.py')
parser.add_argument("-b", "--blast", nargs=1, required=True, help="full path of blast folder", metavar="PATH")
parser.add_argument("-u", "--uniprot", nargs=1, required=True, help="full path of uniprot file", metavar="PATH")
args = parser.parse_args()
print(args)
onlyfiles = [ f for f in os.listdir(args.blast[0]) if os.path.isfile(os.path.join(args.blast[0],f)) ]
print(args.uniprot[0])
for f in onlyfiles:
    
    print(f)
    #print(args.uniprot[0])
    cmd="python D:/Code/Proteomics/Python/IsoformsSAP/UniProteinLocation.py -b "+os.path.join(args.blast[0],f)+" -u "+args.uniprot[0]+" -o "+os.path.join(args.blast[0],"Location",f)
    print(cmd)
    os.system(cmd)
    #break

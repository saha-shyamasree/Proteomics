## This python code runs PITORF:All for a fasta file

import argparse
import subprocess as sp

parser = parser = argparse.ArgumentParser(description='This python code runs PITORF:All for a fasta file')
parser.add_argument("-i", "--input", help="full path of fasta", metavar="String")
parser.add_argument("-o", "--output", help="full path of of output file", metavar="String")
parser.add_argument("--minMain", help="the minimum length of the main ORF, the default=200", metavar="Integer")
parser.add_argument("--start", help="the start codon used to search for ORFs, the default value is M (methionine, corresponding to codon ATG)", metavar="Character")
parser.add_argument("--stop", help="the stop codon used to search for ORFs, the default value is - (corresponding to codons TAA, TAG and TGA)", metavar="String")
parser.add_argument("--uorf5output", help="the file name of the list of 5' uORFs, this also indicates that 5 uORFs is required", metavar="String")
parser.add_argument("--max5", help="the maximum length of 5' uORFs, only works when -uorf5output is set", metavar="Integer")
parser.add_argument("--min5", help="the minimum length of 5' uORFs, only works when -uorf5output is set", metavar="Integer")
parser.add_argument("--distance5", help="the maximum distance between the end position of 5' uORF and the start position of its corresponding main ORF, only works when -uorf5output is set. When this value is set to 0, means that uORF must overlap with the main ORF", metavar="Integer")
parser.add_argument("--uorf3output", help="the file name of the list of 3' uORFs, this also indicates that 3 uORFs is required", metavar="String")
parser.add_argument("--min3", help="the minimum length of 3' uORFs, only works when -uorf3output is set", metavar="Integer")

args = parser.parse_args()
print(args)
#subprocess.call(["perl /data/home/btw796/Code/Proteomics/Perl/orfall.pl", "-h"])
sp.check_output(["perl", "D:\Code\Proteomics\Perl\orfall.pl",args.input,args.output])

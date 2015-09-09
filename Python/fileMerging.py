### This is a python script to merge 2 files. This code is used to merge ORFs file and a contaminants list.

import argparse

parser = parser = argparse.ArgumentParser(description='This is a python script to merge two files.')
parser.add_argument("-1", "--file1", help="full path of file1", metavar="File path")
parser.add_argument("-2", "--file2", help="full path of file2", metavar="File path")
parser.add_argument("-o", "--out", help="full path of of output file", metavar="File path")

args = parser.parse_args()
merge(args.file1, args.file2, args.out)

def merge(file1, file2, out):
    with open(file1,'r') as f1, open(file2,'r') as f2, open(out, 'w') as o:
        for line in f1:
            o.write(line)
        for line in f2:
            o.write(line)
    
    


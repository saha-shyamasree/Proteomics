import sys
import re
import argparse

def find_replace(regex, filename,replace):
    handle = open(filename,"rU")
    
    for line in handle:
        line=line.rstrip()
        if line.startswith(">"):
            print(re.sub(regex,r'\1',line))
        else:
            print(line);
def manipulate(regex, filename, replace):
    handle = open(filename, "rU")
    for line in handle:
        if line.startswith("c"):
            #print(line)
            matches = re.findall(regex, line)
            if len(matches)>0:
                ids=""
                for i in range(len(matches)):
                    ids=ids+";"+matches[i]
                if not ids.startswith(">"):
                    ids=">"+ids.lstrip(';');
                ids=ids.lstrip(';');
                print(ids)
            else:
                line=line.rstrip()
                print(line)
        else:
            line=line.rstrip();
            print(line)

parser = argparse.ArgumentParser(description='Manipulate fasta header.')
parser.add_argument("filename", help="fasta file name.")
parser.add_argument("--regex_index", type =int,  choices=[0, 1, 2, 3], help="\nregex indices are \n0 = User Define\n1=Uniprot\n2=Everything after > (deafault)\n3=Upto first space)")
parser.add_argument("--regex",help="If regex_index is 0, write your own regular expression to extract pattern from fasta header.")
parser.add_argument("--replace",help="If you want custom text to replace regex string put this here.")
args = parser.parse_args()

if args.regex_index == 0:
    if args.regex:
        if args.replace:
            find_replace(args.regex, args.filename,args.replace)
        else:
            manipulate(args.regex, args.filename)
    else:
        raise Exception("No regular expression given for custom pattern")
elif args.regex_index == 1:
    manipulate(">.*\|(.*)\|",args.filename)
elif args.regex_index == 2:
    manipulate(">([^ ]*)", args.filename)
elif args.regex_index == 3:
    manipulate(">([^\|]*)", args.filename)
else:
    manipulate(">([^ ]*)",args.filename)

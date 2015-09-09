## This python code generate reverse decoy sequences.

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-1", "--file1", help="full path of file1", metavar="String")
parser.add_option("-o", "--out", help="full path of of output file", metavar="String")

(options, args) = parser.parse_args()
decoy(options.file1, options.out)

def decoy(file1, out):
    with open(file1,'r') as f1, open(out, 'w') as o:
        prevID=None
        for line in f1:
            if line.startswith(">"):
                o.write(line)
                prevID=line
            else:
                o.write(line)
                if prevID:
                    o.write(line+"_REVERSED")
                    o.write(line[::-1])


import sys
count=0
handle = open(sys.argv[1], "rU")
for line in handle:
    if line.startswith(">"):
        count=count+1
print(count);
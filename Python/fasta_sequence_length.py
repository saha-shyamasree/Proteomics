import sys
import re
from Bio import SeqIO

#seqLength={}
handle = open("human_adenovirus.fasta", "rU")
wrHandle = open("uniprot_human_adeno_length.csv", "w")
for record in SeqIO.parse(handle, "fasta") :
    rid=re.sub(r'^(sp|tr)\|([^\|]+)\|(.*)?',r'\2',record.id)
    #print(str(len(record)))
    wrHandle.write(rid+","+str(len(record))+"\n")
    #seqLength[record.id]=len(record)
    
handle.close()
wrHandle.close()
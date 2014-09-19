from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bs4 import BeautifulSoup
import time
import http.client
import sys
import re
import argparse
import csv

def parseExpasyData(data):
    #pattern=re.compile("(.*)?Number of amino acids:(.*)?Molecular weight:(.*)?Total number of negatively charged residues \(Asp \+ Glu\):(.*)?")
    #pattern=re.compile("Number")
    data=re.sub('\s+',' ',data).strip()
    #print(data)
    match=re.findall("Number of amino acids:(.*)?Molecular weight:(.*)?Theoretical pI:(.*)?Amino acid composition:(?:.*)?The estimated half-life is:(.*)?Instability index: The instability index \(II\) is computed to be(.*)?This classifies the protein as(.*)?Aliphatic index:(.*)?Grand average of hydropathicity \(GRAVY\):(.*)?sib_body",data,re.DOTALL);
    
    if match:
        for i in range(0,len(match)):
            match[i] = [x.strip(' ') for x in match[i]]
            return(match[i])
    else:
        match=re.findall("Number of amino acids:(.*)?Molecular weight:(.*)?Theoretical pI:(.*)?Amino acid composition:(?:.*)?Estimated half-life:(?:.*)?(half-life can not be computed\.) Instability index: The instability index \(II\) is computed to be(.*)?This classifies the protein as(.*)?Aliphatic index:(.*)?Grand average of hydropathicity \(GRAVY\):(.*)?sib_body",data,re.DOTALL);
        if match:
            for i in range(0,len(match)):
                match[i] = [x.strip(' ') for x in match[i]]
                return(match[i])
        else:
            print("match not found")

parser = argparse.ArgumentParser(description='Calculating Protein half-life and other properties')
parser.add_argument("filename", help="fasta file name")
args = parser.parse_args()

handle = open(sys.argv[1], "rU")
print("ProteinID\tNumber of amino acids\tMolecular weight\tTheoretical pI\tHalf-life\tThe instability index\tProtein classification\tAliphatic index\tGrand average of hydropathicity")
if sys.argv[1].endswith('fasta'):
    for record in SeqIO.parse(handle, "fasta"):
        conn = http.client.HTTPConnection("web.expasy.org")
        conn.request("GET", "/cgi-bin/protparam/protparam/?sequence="+record.seq)
        res = conn.getresponse()
        data = res.read()
        #print(data)
        parsedData = BeautifulSoup(data)
        bodyData=parsedData.body.find('div', attrs={'id':'sib_body'}).text
        dataList=parseExpasyData(bodyData)
        dataList.insert(0,record.id)
        print('\t'.join(dataList))
        #print(type(bodyData).__name__)
        #print(bodyData)
        #print(record.id+","+halfLife)
        time.sleep(3)
else:
    c=0
    for line in csv.reader(handle,delimiter='\t'):
        if c>0:
            c=c+1
            if c>0:
                conn = http.client.HTTPConnection("web.expasy.org")
                uId=re.search("\|(.*)?\|",line[2])
                if uId:
                    uId=uId.group(1)
                    uId=re.sub("-\d+","",uId)
                    #print(uId)
                    #print("/cgi-bin/protparam/protparam/?"+uId+"@noft@");
                    #conn.request("GET", "/cgi-bin/protparam/protparam1/?Q9Y2G9@noft@")
                    conn.request("GET", "/cgi-bin/protparam/protparam1/?"+uId+"@noft@")
                    res = conn.getresponse()
                    data = res.read()
                    #print(data)
                    parsedData = BeautifulSoup(data)
                    bodyData=parsedData.body.find('div', attrs={'id':'sib_body'}).text
                    #print(bodyData)
                    dataList=parseExpasyData(bodyData)
                    if dataList:
                        dataList.insert(0,uId)
                        print('\t'.join(dataList))
                        #print(type(bodyData).__name__)
                        #print(bodyData)
                        #print(record.id+","+halfLife)
                        time.sleep(5)
                    else:
                        print(uId+ ": dataList of type none\n\n")
                    #break
                else:
                    print("Uniprot id not found")
                conn.close()
        else:
            c=1
        #if c==2:
        #    break

handle.close
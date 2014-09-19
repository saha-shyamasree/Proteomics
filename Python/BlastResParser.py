#!/usr/bin/python
# -*- coding: utf-8 -*-

import csv
import sys
import re
import glob
import os

fName=r"C:\Users\shyama\Dropbox\results\trinity_uniprot_blast.csv"
reader = open(fName, 'rU')
for line in reader:
    line=line.rstrip();
    if line.startswith('comp'):
        matches = re.findall("([^,]+)", line)
        if len(matches)>0:
            #print(matches[0])
            ids=""
            t_ids=re.findall("comp[^ ]+",matches[0])
            if len(t_ids)>0:
                for i in range(len(t_ids)):
                    ids=ids+";"+t_ids[i]
                ids=ids.lstrip(';')
                print(line.replace(matches[0],ids))
    else:
        print(line)
                
            
#te=BlastResParser(fName)
#te.read()
## This code reads pasa_run.log.dir/alignment_assembly_subclustering.out from PASA output folder and creates a tab delimated file with 2 columns, 1. cluster number and 2. assemblies of the cluster. A cluster represents overlapping assemblies from PASA.
import re

def readFile(filePath, wrtPath):
    clusterIdStr="// Processing cluster: ";
    clusterAlignments="Individual Alignments: (";
    subCluster="sub-cluster: "
    clusterString="";
    
    with open(filePath, 'r') as RD, open(wrtPath, 'w') as WRT:
        WRT.write("Cluster Id, Alignment Count, Assemblies\n")
        flag=0
        count=0
        for line in RD:
            if line.startswith(clusterIdStr):
                clusterId=re.findall('\d+',line)
                flag=1;
                if clusterId!=None:
                    if count==0:
                        WRT.write(clusterId[0]+",")
                    else:
                        WRT.write("\n"+clusterId[0]+",")
                count=count+1
            elif line.startswith(clusterAlignments):
                alignmntsCount=re.findall('\d+',line)
                if alignmntsCount!=None:
                    WRT.write(alignmntsCount[0]+",")
            elif line.startswith(subCluster):
                line=line.strip();
                assemblies=line[line.index("asmbl_"):]
                if flag==1:
                    WRT.write(assemblies);
                else:
                    WRT.write("|"+assemblies);
                flag=0

clusterInfoFile="D:/data/PASA/human_adeno/Output/alignment_assembly_subclustering.out"
clusterCSV="D:/data/PASA/human_adeno/Output/PostProcessing/PASAAssemblyCluster.csv"

readFile(clusterInfoFile, clusterCSV)
    
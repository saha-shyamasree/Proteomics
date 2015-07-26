## this is the python code, re-implementation of my R code, to compute peptide and protein overlaps, comparing peptide/protein
## identification results from MSGF+ (identification) and mzIdenML-lib(post-processing)

import pandas as pd;

def readDelimetedFile(filename, sep):
    dataFrameObj = pd.read_table(filename, sep=sep)
    return dataFrameObj;

#Mat is a DataFrame object containing protein grouping results from MSGF+
def ProteinGroupFiltered(Mat, rev, peptide, pepThreshold):
    ##removing proteins from decoy DB.
    if rev==1:
        

def AnchorProtein(Mat):
    
    
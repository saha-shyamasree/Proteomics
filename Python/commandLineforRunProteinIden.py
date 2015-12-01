#### Calling runProteinIdentification.py on Desktop
import os
msPath="D:/data/Oliver/"
transPath="D:/data/Oliver/PASA/"
transToMSMap={
#    'G11':'E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S1.RAW.-1.mgf',
#    'G15':'E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S2.RAW.-1.mgf',
#    'G17':'E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S3.RAW.-1.mgf',
#    'G29a':'E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S4.RAW.-1.mgf',
#    'G30':'E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S5.RAW.-1.mgf',
#    'G75':'E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S19.RAW.-1.mgf',
#    'G102':'E_Vinni_Vinni_projects_b1368p25_b1368p25_ECM_S24.RAW.-1.mgf',
#    'G138':'3eab50639f039213___b1368p92_ECM_S30.RAW.-1.mgf'
}

for key,value in transToMSMap.items():
    os.system("python runProteinIdentificationAndPostProcessing.py -i "+msPath+value+
              " -o "+transPath+key+"/"+key+".assemblies.fasta.transdecoder.pep.mzid -d "+
              transPath+key+"/"+key+".assemblies.fasta.transdecoder.pep -m"+
              msPath+"modifications.txt")
    

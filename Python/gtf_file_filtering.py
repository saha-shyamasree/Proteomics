import sys;
import re;
if len(sys.argv)==2:
    with open(sys.argv[1]) as f:
        proteinCodingStr = ["\"protein_coding\";","\"polymorphic_pseudogene\";","\"IG_C_gene\";","\"IG_D_gene\";","\"IG_J_gene\";","\"IG_V_gene\";","\"TR_C_gene\";","\"TR_D_gene\";","\"TR_J_gene\";","\"TR_V_gene\";"]
        for line in f:
            rLine=re.sub("\t"," ",line.rstrip('\n'));
            fields=rLine.split(" ");
            if len(fields)>=3 and fields[2]=="gene":
                gidx=fields.index('gene_biotype')
                ##gidx is index of gene_biotype, hence gidx+1 is index of the value. 
                if len(fields)>=gidx+2 and fields[gidx+1] in proteinCodingStr:
                    part1="\t".join(fields[0:8])
                    part2=" ".join(fields[8:])
                    all="\t".join([part1,part2])
                    print(all)
            else:
                try:
                    ##tidx is index of transcript_biotype, hence tidx+1 is index of the value and size of the fields variable should be tidx+2
                    tidx=fields.index('transcript_biotype')
                    print(tidx)
                    if len(fields)>=tidx+2 and fields[tidx+1] in proteinCodingStr:
                        part1="\t".join(fields[0:8])
                        part2=" ".join(fields[8:])
                        all="\t".join([part1,part2])
                        print(all)
                except ValueError:
                    continue;
                
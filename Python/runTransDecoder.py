##This code runs transDecoder
import argparse
import os

parser = argparse.ArgumentParser(description='Fasta file filtering based on a header list given')
parser.add_argument("-g","--gff3", nargs=1, required=True, help="gff3 file name")
parser.add_argument("-t","--transcripts", nargs=1, required=True, help="Transcript fasta file name")
parser.add_argument("-p","--transDPath", default=['/data/home/btw796/Prog/PASApipeline-master/scripts/'], help="Transdecoder path")
parser.add_argument("-o","--out", nargs=1, help="Output filename")

args = parser.parse_args()
print(args)
print(args.out)
if args.out==None:
    args.out=args.transcripts
#print(args.transcripts)

command=" ".join(["echo",args.transDPath[0], "pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta", args.transcripts[0],"--pasa_transcripts_gff3", args.gff3[0], "\" | qsub -cwd -V -l h_vmem=8G -l h_rt=24:0:0"])
print(command)
print(args.transcripts)

#echo "/data/home/btw796/Prog/PASApipeline-master/scripts/pasa_asmbls_to_training_set.dbi -m 11 --pasa_transcripts_fasta $sample.assemblies.fasta --pasa_transcripts_gff3 $sample.pasa_assemblies.gff3" | qsub -cwd -V -l h_vmem=8G -l h_rt=24:0:0
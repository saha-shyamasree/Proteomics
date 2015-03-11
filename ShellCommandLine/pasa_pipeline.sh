################ Cleaning transcript sequences ########################
./seqclean ../../human_adeno_data/human_adenovirus_trinity_assembled_transcripts.fasta -o ../../human_adeno_data/

################ Filtering Protein coding genes out from the Ensembl GTF file #######################
python ~/Code/Proteomics/Python/gtf_file_filtering.py Homo_sapiens.GRCh38.78.gtf > Homo_sapiens.GRCh38.78_Protein_Coding.gtf
################ Converting GTF file to GFF3 Format, using gtf2gff tool from Apocrita Cluster ###############
gtf2gff.pl < Homo_sapiens.GRCh38.78_Protein_Coding.gtf --gff3 --out=Homo_sapiens.GRCh38.78_Protein_Coding.gff3

################ Validating the GFF3 file for PASA ######################
../misc_utilities/pasa_gff3_validator.pl Homo_sapiens.GRCh38.78_Protein_Coding.gff3

############### Run the PASA alignment assembly pipeline like so: ###################

../scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g genome_sample.fasta -t human_adenovirus_trinity_assembled_transcripts_clean.fasta -T -u all_transcripts.fasta --TRANSDECODER --ALIGNERS blat,gmap --CPU 2 --MAX_INTRON_LENGTH 300000 --ALT_SPLICE --stringent_alignment_overlap 30.0 --gene_overlap 50.0





-c configuration_file [DONE]
--ALIGNERS 'gmap,blat' [DONE]
--MAX_INTRON_LENGTH 300000 [DONE]
--cufflinks_gtf <filename>      :incorporate cufflinks-generated transcripts
-C #flag, create MYSQL database [DONE]
-R [DONE]
-A
--ALT_SPLICE [DONE]
-g <filename> #  genome sequence FASTA file (should contain annot db asmbl_id as header accession.) [DONE]
-t <transcript_file> [DONE]
--TDN <filename> file containing a list of accessions corresponding to Trinity (full) de novo assemblies (not genome-guided)
-T               #flag,transcript db were trimmed using the TGI seqclean tool. [DONE]
-u <filename>   #value, transcript db containing untrimmed sequences (input to seqclean) [DONE]
--TRANSDECODER  #flag, run transdecoder to identify candidate full-length coding transcripts [DONE]
--stringent_alignment_overlap 30.0  (#suggested: 30.0)  overlapping transcripts must have this min % overlap to be clustered.
--gene_overlap 50.0  #(suggested: 50.0)  transcripts overlapping existing gene annotations are clustered.  Intergenic alignments are clustered by default mechanism.
-L
--annots_gff3 <filename>  #existing gene annotations in gff3 format.




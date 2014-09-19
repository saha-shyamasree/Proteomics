##### remember to remove all commas from fasta header before following steps. ###############

makeblastdb -in C:\Users\shyama\Dropbox\results\blastdb\uniprot_P_alecto_hendra.fasta -dbtype prot -out C:\Users\shyama\Dropbox\results\blastdb\uniprot_P_alecto_hendra

blastp -db C:\Users\shyama\Dropbox\results\blastdb\uniprot_P_alecto_hendra -query C:\Users\shyama\Dropbox\results\blastdb\Trinity_P_alecto-ORF_header2.fasta -gapopen 6 -gapextend 2 -out C:\Users\shyama\Dropbox\results\blastdb\uni_db_blast -outfmt 5
blastp -db C:\Users\shyama\Dropbox\results\blastdb\Trinity_P_alecto-ORF_header2 -query C:\Users\shyama\Dropbox\results\blastdb\uniprot_P_alecto_hendra.fasta -gapopen 6 -gapextend 2 -out C:\Users\shyama\Dropbox\results\blastdb\p_alecto_blast -outfmt 5

makeblastdb -in C:\Users\shyama\Dropbox\results\blastdb\Trinity_P_alecto-ORF_header.fasta -dbtype prot -out C:\Users\shyama\Dropbox\results\blastdb\Trinity_P_alecto-ORF_header
blastp -db C:\Users\shyama\Dropbox\results\blastdb\Trinity_P_alecto-ORF_header2 -query C:\Users\shyama\Dropbox\results\blastdb\uniprot_P_alecto_hendra.fasta -gapopen 6 -gapextend 2 -out C:\Users\shyama\Dropbox\results\blastdb\p_alecto_blast -outfmt 5

blastp -db E:\Data\HUMAN\database\human_adenovirus.fasta -query E:\Data\HUMAN\database\Trinity\Supplementary_dataset_1_-_human_trinity_assembled_transcripts-ORF.fasta -gapopen 6 -gapextend 2 -out G:\trinity_PITORF_blast -outfmt 5


makeblastdb -in C:\Users\Shyamasree\Dropbox\results\Bat_human_hendra\Human\Trinity_HeV293-ORF.fasta -dbtype prot -out C:\Users\Shyamasree\Dropbox\results\blastdb\trinity_human_hendra
blastp -db C:\Users\Shyamasree\Dropbox\results\blastdb\trinity_human_hendra -query C:\Users\Shyamasree\Dropbox\results\Bat_human_hendra\Human\uniprot.fasta -gapopen 6 -gapextend 2 -out I:\Human-Hendra\trinity_db_blast -outfmt 5

makeblastdb -in C:\Users\Shyamasree\Dropbox\results\Bat_human_hendra\Human\uniprot.fasta -dbtype prot -out C:\Users\Shyamasree\Dropbox\results\blastdb\uniprot_human_hendra
blastp -db C:\Users\Shyamasree\Dropbox\results\blastdb\uniprot_human_hendra -query C:\Users\Shyamasree\Dropbox\results\Bat_human_hendra\Human\Trinity_HeV293-ORF.fasta -gapopen 6 -gapextend 2 -out I:\Human-Hendra\uniprot_db_blast -outfmt 5
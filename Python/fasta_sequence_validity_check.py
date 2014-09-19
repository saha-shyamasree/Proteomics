with open("E:\Data\HUMAN\database\Trinity\human_trinity_ORFS_windows.fasta", "wt") as fout:
    with open("E:\Data\HUMAN\database\Trinity\Supplementary dataset 2 - human trinity ORFS.fasta", "rt") as fin:
        for line in fin:
            fout.write(line.replace('\r\n', '\n'))
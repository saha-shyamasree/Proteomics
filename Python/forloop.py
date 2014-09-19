import csv

with open(r"C:\Users\shyama\Dropbox\Gideon\sORF_list_for_python.csv", "r") as input_file_1, \
    open(r"C:\Users\shyama\Dropbox\Gideon\test.csv", "r") as input_file_2:

    sORF_file = list(csv.reader(input_file_1, delimiter = ";", dialect = "excel"))
    test_file = list(csv.reader(input_file_2, delimiter = ";", dialect = "excel"))
    #output_file = csv.writer(output, delimiter = ",", dialect = "excel")
    count = 0
    prev = 0
    for sORF_row in sORF_file: #csv.reader(input_file_1, delimiter = ";", dialect = "excel"):
        count += 1
        temp_sORF_row = sORF_row
        #print(temp_sORF_row)
        for evalue_row in test_file: #csv.reader(input_file_2, delimiter = ";", dialect = "excel"):
            #results_row = sORF_row
            print(sORF_row[1], evalue_row[0])
            if evalue_row[0] in sORF_row[1]:
                print("it works")
                 #results_row.append(evalue_row[0])
                 #results_row.append(evalue_row[2])
                 #results_row.append(evalue_row[3])
                 #output_file.writerow(results_row)
            else:
                if count == prev:
                    print("it's not working", count)
                    prev = count
        #print(sORF_row)
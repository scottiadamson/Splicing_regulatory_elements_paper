import csv

b = open('/home/CAM/adamson/Recount2_files/sample_ids.tsv', 'w')
writer = csv.writer(b, delimiter = '\t')
i = 0
for file_type in ['SRP012682', 'TCGA']:
    a = open('/home/CAM/adamson/Recount2_files/' + file_type + '/' + file_type +'.tsv', 'r')
    reader = csv.DictReader(a, delimiter = '\t')
    for line in reader:
        writer.writerow([i, file_type, line['run']])
        i+=1
    a.close()
b.close()


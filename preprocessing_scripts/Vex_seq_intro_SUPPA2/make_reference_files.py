import csv, sys

base_dir = '/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test/'
a = open(base_dir + 'ref2var.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
var2ref = {}
for line in reader:
    for variant in line[1].split(','):
        var2ref[variant] = line[0]
    var2ref[line[0]] = line[0]
a.close()

b = open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/all.fa', 'r')
seqs = {}
for line in b:
    line = line.replace('\n', '')
    if line[0] == '>':
        seq_id = line.replace('>', '')
    else:
        seqs[seq_id] = line
b.close()

variant = sys.argv[1]

c = open(base_dir + 'references/' + variant + '.fa', 'w')
c.write('>' + variant + '\n' + seqs[variant][2:] + '\n')
c.close()

#10000  .   exon    1   46  .   +   .   transcript_id "10000_excluded"; gene_id "10000";
#10000   .   exon    1005    1187    .   +   .   transcript_id "10000_excluded"; gene_id "10000";
#10000   .   exon    1   46  .   +   .   transcript_id "10000_included"; gene_id "10000";
#10000   .   exon    241 285 .   +   .   transcript_id "10000_included"; gene_id "10000";
e = open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/all.gtf', 'r')
reader = csv.reader(e, delimiter = '\t')
f = open(base_dir + 'references/' + variant + '.gtf', 'w')
writer = csv.writer(f, delimiter = '\t',quotechar='', escapechar= '\\', quoting=csv.QUOTE_NONE)
i=0; these_lines = {} 
for line in reader: 
    if line[0] == variant:
        these_lines[i] = line
        i+=1
max_gene = int(these_lines[1][4])-2
writer.writerow([variant, '.', 'gene', 1, max_gene, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant+ '"; '])])
writer.writerow([variant, '.', 'transcript', 1, max_gene, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant + '"', ' transcript_id "' + variant + '_excluded";'])])
writer.writerow([variant, '.', 'exon', 1 , int(these_lines[0][4]) -2, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant+ '"', ' transcript_id "' + variant + '_excluded";'])])
writer.writerow([variant, '.', 'exon', int(these_lines[1][3]) - 2, int(these_lines[1][4]) - 2, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant+ '"', ' transcript_id "' + variant + '_excluded";'])]) 
writer.writerow([variant, '.', 'transcript', 1, max_gene, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant+ '"', ' transcript_id "' + variant + '_included";'])])
writer.writerow([variant, '.', 'exon', 1 , int(these_lines[2][4]) - 2, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant+ '"', ' transcript_id "' + variant + '_included";'])])
writer.writerow([variant, '.', 'exon', int(these_lines[3][3]) - 2, int(these_lines[3][4]) - 2, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant+ '"', ' transcript_id "' + variant + '_included";'])]) 
writer.writerow([variant, '.', 'exon', int(these_lines[4][3]) - 2, int(these_lines[4][4]) - 2, '.', '+', '.', ';'.join(['gene_id ' + '"' + variant+ '"', ' transcript_id "' + variant + '_included";'])]) 
e.close();f.close()


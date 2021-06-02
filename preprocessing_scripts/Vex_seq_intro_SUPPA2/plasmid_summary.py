import csv
from collections import defaultdict

base_dir = '/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/'
barcode_ids = set()
a = open(base_dir + 'barcode_samples.txt', 'r')
for line in a:
    barcode_ids.add(line.replace('\n', ''))
a.close()

mapped_depth = {'1o': defaultdict(lambda:0), '2o': defaultdict(lambda:0)}
for sample in ['1o', '2o']:
    a = open(base_dir + sample + '_uniquely_assigned.txt', 'r')
    for line in a:
        barcode_id = line.split(' ')[0]
        mapped_depth[sample][barcode_id] = line.strip().split(' ')[1]
    a.close() 

def calc_emp_AF(line):
    max_AF = 0
    DP = int(line[9].split(':')[2])
    for AD in line[9].split(':')[3].split(',')[1:]:
        max_AF = max(max_AF, float(AD)/DP)
    return max_AF

variant_max_AF = {'1o': {}, '2o': {}}
for sample in ['1o']:
    #variants identified
    b=open(base_dir + 'vcfs/' + sample + '_merged.vcf2', 'r')
    ##CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  /home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/alignments/1o_SRE/7000_AACTGGCCGA_sorted.bam
    #10002_TACTCGGGGA        110     .       G       T       32.1827 .       DP=58;VDB=0.0221621;SGB=-0.511536;RPB=0.994247;MQB=1;BQB=0.911825;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=26,0,3,0;MQ=60     GT:PL:DP:AD     0/1:69,0,255:29:26,3
    reader = csv.reader(b, delimiter = '\t')
    for line in reader:
        if line[0][0] != '#':
            if line[0] not in variant_max_AF:
                variant_max_AF[sample][line[0]] = calc_emp_AF(line)
            else:
                variant_max_AF[sample][line[0]] = max(variant_max_AF[sample][line[0]], calc_emp_AF(line))
    b.close()

b = open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/plasmid_statistics2.tsv', 'w')
writer = csv.writer(b, delimiter = '\t')
writer.writerow(['barcode', '1o_depth', '1o_mapped', '1o_variant_max_AF', '2o_depth', '2o_mapped', '2o_variant_max_AF'])
#writer.writerow(['barcode', '1o_mapped', '1o_max_variant_count', '1o_max_variant_AF', '2o_mapped', '2o_variant_max_AF'])
for barcode in sorted(barcode_ids):
    line = [barcode]
    line.append(sum(1 for line in open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/sorted_fastqs/1o_SRE/' + barcode + '_R1.fastq'))/4)
    line.append(mapped_depth['1o'][barcode])
    if barcode in variant_max_AF['1o']:
        line.append(variant_max_AF['1o'][barcode])
    else:
        line.append(0)
    line.append(sum(1 for line in open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/sorted_fastqs/2o_SRE/' + barcode + '_R1.fastq'))/4)
    line.append(mapped_depth['2o'][barcode])
    if barcode in variant_max_AF['2o']:
        line.append(variant_max_AF['2o'][barcode])
    else:
        line.append(0)
    writer.writerow(line)
b.close()




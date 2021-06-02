import csv, gzip, sys
from datetime import datetime
from collections import defaultdict
print(datetime.now())

#note that this is neccessary so that this pipeline can be run in pieces, otherwise too many files are open at once and it makes computers (and your HPC admin) very sad. 
a = open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/SRE_pool_sequence_no_dups.tsv', 'r')
reader = csv.DictReader(a,delimiter = '\t')
#group_id   test_sequence   full_sequence   barcode
#1   CAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTG  CTGACTCTCTCTGCCTCCAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTGCAATTGACTACTAGTGTCCCCTATCTCTAGAGGGCCCGTTTA   GTCCCCTATC
BC2id = {}
for line in reader:
        BC2id[line['barcode']] = line['group_id']
a.close()



b = open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/plasmid_statistics2.tsv', 'r')
reader = csv.reader(b, delimiter = '\t')
#barcode    1o_depth    1_mapped    1o_variant_max_AF   2o_depth    2_mapped    2o_variant_max_AF
#10000_AGGGGCTAAG    14  32  0   6   12  0
bad_BCs = set()
next(reader)
for line in reader:
    if float(line[3]) > 0.2:
        bad_BCs.add(line[0].split('_')[1])
b.close()

def rev_comp(seq):
    comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' :'N'}
    return ''.join([comp[x] for x in seq])[::-1]

sample = sys.argv[1]
ref_base = sys.argv[2]

#0-100_R1.fastq.gz
#0-100_R2.fastq.gz
lower_range = int(ref_base.split('-')[0])
upper_range = int(ref_base.split('-')[1])

total_reads = 0; assigned_reads = 0; bad_barcodes = 0
files = {};Read_name_group = {}
for seq in range(lower_range, upper_range):
    files[str(seq)] = gzip.open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/temp_fastqs/' + sample + '/' + str(seq) + '_R2.fastq.gz', 'wt')
a = gzip.open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/sorted_fastqs/' + sample + '/' + ref_base + '_R2.fastq.gz', 'rt')
i=0
for line in a:
    if i == 0:
        name = line
    if i == 1:
        seq = line
    if i == 2:
        strand = line
    if i == 3:
        qual = line
        RC_seq = rev_comp(seq[0:-1])
        barcode = RC_seq[RC_seq.find('CTAAGCTCGCTCTAGT') + 16: RC_seq.find('CTAAGCTCGCTCTAGT') + 26]
        total_reads += 1
        group_id = BC2id[barcode]
        if barcode not in bad_BCs:
            files[group_id].write(name + seq + strand + qual)
            Read_name_group[name.split(' ')[0]] = group_id
        else:
            bad_barcodes +=1
        i=-1
    i+=1
a.close()
print('total:    ' + str(total_reads))
#print('assigned: ' + str(assigned_reads))
#print('bad BCs:  ' + str(bad_barcodes))
print(datetime.now())
for seq in range(lower_range, upper_range):
    files[str(seq)].close()
    files[str(seq)] = gzip.open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/temp_fastqs/' + sample + '/' + str(seq) + '_R1.fastq.gz', 'wt')
#R1
a = gzip.open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/sorted_fastqs/' + sample + '/' + ref_base + '_R1.fastq.gz', 'rt')
i=0
for line in a:
    if i == 0:
        name = line
    if i == 1:
        seq = line
    if i == 2:
        strand = line
    if i == 3:
        qual = line
        if name.split(' ')[0] in Read_name_group:
            group_id = Read_name_group[name.split(' ')[0]]
            files[group_id].write(name + seq + strand + qual)
        i=-1
    i+=1
a.close()
for seq in range(lower_range, upper_range):
    files[str(seq)].close()

print(datetime.now())

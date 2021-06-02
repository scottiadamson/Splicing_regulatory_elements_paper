import csv, gzip, sys
from datetime import datetime
print(datetime.now())
print(sys.argv)
#group_id    test_sequence   full_sequence   barcode
#1   CAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTG  CTGACTCTCTCTGCCTCCAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTGCAATTGACTACTAGTGTCCCCTATCTCTAGAGGGCCCGTTTA   GTCCCCTATC
ref2var = {}
a = open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test/ref2var.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
for line in reader:
    ref2var[line[0]] = []
    for var in line[1].split(','):
        ref2var[line[0]].append(var)
    ref2var[line[0]].append(line[0])
a.close()

sample = sys.argv[1]
reference_group = sys.argv[2]

clip_1 = set(['K562_2', 'HepG2_2', 'K562_NT_2', 'K562_UPF1_1'])
clip_2 = set(['K562_3', 'HepG2_3', 'K562_UPF1_2'])
def clip_read(sample, seq, qual):
    if sample in clip_1:
        seq = seq[1:]
        qual = qual[1:]
    elif sample in clip_2:
        seq = seq[2:]
        qual = qual[2:]
    return seq, qual

for read_type in ['R1', 'R2']:
    files = {}
    for variant in ref2var[reference_group]:
        files[variant] = gzip.open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test/temp_fastqs/' + sample + '/' + variant + '_' + read_type + '.fastq.gz', 'wt')
    a = gzip.open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test/sorted_fastqs/' + sample + '/' + reference_group + '_' + read_type + '.fastq.gz', 'rt')
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
            group_id = name.split(':')[-1].replace('\n', '')
            if read_type == 'R1':
                seq, qual = clip_read(sample, seq, qual)
            files[group_id].write(name + seq + strand + qual)
            i=-1
        i+=1
    a.close()
    for variant in ref2var[reference_group]:
        files[variant].close() 

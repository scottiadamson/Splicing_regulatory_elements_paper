import csv
from collections import defaultdict

samples = ['HepG2', 'K562'] #, 'K562_NT', 'K562_UPF1']
base_dir = '/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/NextSeq_Sept_2019/transcriptome_assembled/'

#1014	1014	44_1003	1014_excluded	216,228,144	216,228,144
junctions = defaultdict(set)
total_counts = {}
for sample in samples:
    reference_seen = set()
    replicates = [sample + '_' + rep for rep in ['1', '2', '3']]
    for replicate in replicates:
        total_counts[replicate] = 0
    a = open(base_dir + sample + '_all_counts.tsv', 'r')
    reader = csv.reader(a, delimiter = '\t')
    for line in reader:
        reference_junction = line[1] + ':' + line[2]
        if reference_junction not in reference_seen:
            for i in range(0,3):
                total_counts[sample + '_' + str(i+1)] += int(line[4].split(',')[i])
            reference_seen.add(reference_junction)
        for i in range(0,3):
                total_counts[sample + '_' + str(i+1)] += int(line[5].split(',')[i])
        if 'novel' in line[3]:
            junctions[line[0]].add(line[2])
    a.close()
print(total_counts)
#10000	.	exon	1003	1185	.	+	.	transcript_id "10000_included"; gene_id "10000";
b = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/NextSeq_Sept_2019/all.gtf', 'r')
reader = csv.reader(b, delimiter = '\t')
c = open(base_dir + 'all_novel_merged.gtf', 'w')
writer = csv.writer(c, delimiter = '\t',quotechar='', escapechar= '\\', quoting=csv.QUOTE_NONE)
seen_exons = defaultdict(set)
junction_rename = {};seen_included = defaultdict(lambda:0)
for line in reader:
    if line[3] != '1':
        line[3] = str(int(line[3]) - 2)
    line[4] = str(int(line[4]) - 2)
    if 'included' in line[8]:
        if seen_included[line[0]] == 0:
            junctions_inc = line[4]
        elif seen_included[line[0]] == 1:
            junctions_inc += '_' + line[3] + '_' + line[4]
        else:
            junctions_inc += '_' + line[3]
        seen_included[line[0]] += 1
        
    writer.writerow(line)

    seen_exons[line[0]].add(line[3] + '_' + line[4])
    #10019	24976	44_219_262_283_290_1003	10019_novel1	0,0,0	3,1,1
    if len(seen_exons[line[0]]) == 3 and int(line[4]) > 1000:
        junctions_exc = junctions_inc.split('_')[0] + '_' + junctions_inc.split('_')[3]
        i =0
        for junction in junctions[line[0]]:
            reformat_junction = ''
            transcript_id = line[0] + '_novel' + str(i)
            junction_rename[line[0] + ':' + junction] = transcript_id
            exon_count = int(len(junction.split('_'))/2 + 1)
            for exon in range(0, exon_count):
                if exon == 0:
                    writer.writerow([line[0], '.', 'exon', 1, junction.split('_')[0], '.', '+', '.', 'transcript_id "' + transcript_id + '"; gene_id "' + line[0] + '";'])
                    reformat_junction += junction.split('_')[0] + '_'
                elif exon != exon_count - 1:
                    writer.writerow([line[0], '.', 'exon', int(junction.split('_')[(exon*2)-1])+1, junction.split('_')[exon*2], '.', '+', '.', 'transcript_id "' + transcript_id + '"; gene_id "' + line[0] + '";'])
                    reformat_junction += str(int(junction.split('_')[(exon*2)-1])+1) + '_' + str(junction.split('_')[exon*2]) + '_'
                else:
                    writer.writerow([line[0], '.', 'exon', int(junction.split('_')[-1]) + 1, line[4], '.', '+', '.', 'transcript_id "' + transcript_id + '"; gene_id "' + line[0] + '";'])
                    reformat_junction += str(int(junction.split('_')[-1]) + 1)
            i += 1
            #if len(junction.split('_')) == 2:
            #    category = classify_event(junctions_exc, reformat_junction)
            #elif len(junction.split('_')) == 4:
            #    category = classify_event(junctions_inc, reformat_junction)
            #else: # these all have 6 splice sites in a 4 exon setup
            #    category = 'additional_exon'
b.close();c.close()


def tpm_transform(line, index, ref): #we don't need to do a length normalization because of the fixed endpoints of the amplicons
    replicate = sample + '_' + str(index+1) #inheriting sample from global in for loop
    if ref:
        count = line[4].split(',')[index] 
    else:
        count = line[5].split(',')[index]
    return round((float(count) *1000000 / total_counts[replicate]), 3)

#1014	1014	44_1003	1014_excluded	216,228,144	216,228,144
#for the output the reference sample TPM will have redundant genes; i.e. the variant and the reference are the same measurement, so the reference samples will have way more counts; hopefully that's okay
for sample in samples:
    a = open(base_dir + sample + '_all_counts.tsv', 'r')
    reader = csv.reader(a, delimiter = '\t')
    b = open(base_dir + sample + '_all_counts_suppa_ref.tsv', 'w')
    ref_writer = csv.writer(b, delimiter = '\t')
    ref_writer.writerow([sample + '_1_ref', sample + '_2_ref', sample + '_3_ref'])
    c = open(base_dir + sample + '_all_counts_suppa_var.tsv', 'w')
    var_writer = csv.writer(c, delimiter = '\t')
    var_writer.writerow([sample + '_1_var', sample + '_2_var', sample + '_3_var'])
    for line in reader:
        if 'novel' in line[3]:
            transcript_id = junction_rename[line[0] + ':' + line[2]]
        else:
            transcript_id = line[3]
        ref_writer.writerow([transcript_id, tpm_transform(line, 0, True), tpm_transform(line, 1, True), tpm_transform(line, 2, True)]) 
        var_writer.writerow([transcript_id, tpm_transform(line, 0, False), tpm_transform(line, 1, False), tpm_transform(line, 2, False)]) 
    a.close();b.close()


import csv
import pandas as pd

sequence_design_info = pd.read_csv('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/take2/SRE_pool_sequences_dups_hg19.tsv', sep = '\t')
#chromosome	position	reference	alternative	intron_lower	intron_upper	exon_start	exon_end	strand	control	test_sequence	full_sequence	barcode	group_id
#chr9	117012834	CG	GT	117012821	117012991	117012894	117012948	+	SRSF1_HepG2|8877_SE	CAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTG	CTGACTCTCTCTGCCTCCAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTGCAATTGACTACTAGTGTCCCCTATCTCTAGAGGGCCCGTTTA	GTCCCCTATC	1
sequence_design_info = dict(zip(sequence_design_info.loc[:,'group_id'], sequence_design_info.to_numpy()))

genotypes = ['HepG2', 'K562'] #, 'K562_NT', 'K562_UPF1']
event_types = ['SE', 'A3', 'A5']

base_dir = '/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/NextSeq_Sept_2019/transcriptome_assembled/'

var2ref = {}
a = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/NextSeq_Sept_2019/ref2var.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
for line in reader:
    for variant in line[1].split(','):
        var2ref[variant] = line[0]
    var2ref[line[0]] = line[0]
a.close()

reference_info = {}
for event_type in event_types:
    #seqname	gene_id	event_id	alternative_transcripts	total_transcripts
    #10000	10000	10000;SE:10000:44-239:283-1003:+	10000_included	10000_included,10000_excluded
    a = open(base_dir + 'master_novel_merged_' + event_type + '_strict.ioe', 'r')
    reader = csv.DictReader(a, delimiter = '\t')
    for line in reader:
        reference_info[line['event_id']] = line
    a.close()

def average_psi(PSIs):
    this_psi = 0
    for psi in PSIs:
        if psi != 'nan' and this_psi != 'nan':
            this_psi += float(psi)
        else:
            this_psi = 'nan'
    if this_psi != 'nan':
        this_psi /= 3
    return this_psi

def add_genomic_locations(row):
    group_id = int(row[0])
    #chr20	57473991	AA	TT	57473911	57474081	57473995	57474040	+	HNRNPL_K562|70862_SE	
    #10000;SE:10000:44-239:283-1003:+
    event_id = row[2].split(':')
    plasmid_junction = False
    if row[3] == 'SE':
        if sequence_design_info[group_id][8] == '+':
            genomic_exon_start = int(event_id[2].split('-')[1]) + sequence_design_info[group_id][4] - 155
            genomic_exon_end = int(event_id[3].split('-')[0]) + sequence_design_info[group_id][4] - 154
        else:
            genomic_exon_start = sequence_design_info[group_id][5] - int(event_id[3].split('-')[0])  + 154
            genomic_exon_end = sequence_design_info[group_id][5] - int(event_id[2].split('-')[1]) + 155
            if int(event_id[3].split('-')[0]) < 154 or int(event_id[2].split('-')[1]) > (155 + 171):
                plasmid_junction = True
            if int(event_id[2].split('-')[0]) != 44 or int(event_id[3].split('-')[1]) not in set([1003, 1004]):
                plasmid_junction = True
        return [sequence_design_info[group_id][0], genomic_exon_start, genomic_exon_end, sequence_design_info[group_id][8], plasmid_junction]
    elif row[3] == 'A3': #10003;A3:10003:44-219:44-239:+
        if sequence_design_info[group_id][8] == '+':
            genomic_long_exon_SS = int(event_id[2].split('-')[1]) + sequence_design_info[group_id][4] - 155
            genomic_short_exon_SS = int(event_id[3].split('-')[1]) + sequence_design_info[group_id][4] - 155
        else:
            genomic_long_exon_SS = sequence_design_info[group_id][5] - int(event_id[3].split('-')[1])  + 155
            genomic_short_exon_SS = sequence_design_info[group_id][5] - int(event_id[2].split('-')[1]) + 155
        if int(event_id[3].split('-')[1]) < 155 or int(event_id[2].split('-')[0]) != 44 or int(event_id[3].split('-')[0]) != 44:
            plasmid_junction = True
    else: #10001;A5:10001:312-1004:262-1004:+
        if sequence_design_info[group_id][8] == '+':
            genomic_long_exon_SS = int(event_id[2].split('-')[0]) + sequence_design_info[group_id][4] - 154
            genomic_short_exon_SS = int(event_id[3].split('-')[0]) + sequence_design_info[group_id][4] - 154
        else:
            genomic_short_exon_SS = sequence_design_info[group_id][5] - int(event_id[3].split('-')[0])  + 154
            genomic_long_exon_SS = sequence_design_info[group_id][5] - int(event_id[2].split('-')[0]) + 154
        if int(event_id[2].split('-')[1]) not in set([1003, 1004]) or int(event_id[3].split('-')[1]) not in set([1003, 1004]) or int(event_id[2].split('-')[0]) > (155 + 171):
            plasmid_junction = True
    return [sequence_design_info[group_id][0], genomic_short_exon_SS, genomic_long_exon_SS, sequence_design_info[group_id][8], plasmid_junction]

#HepG2_ref_A3_1  HepG2_ref_A3_2  HepG2_ref_A3_3  HepG2_var_A3_1  HepG2_var_A3_2  HepG2_var_A3_3
#10003;A3:10003:44-219:44-239:+  0.9935923469157528      1.0     0.990912963225898       0.8115690675372373      0.8301996285979573      0.8332971800433839
for genotype in genotypes:
    c = open(base_dir + genotype + '_suppa_postprocessed.tsv', 'w')
    writer = csv.writer(c, delimiter = '\t')
    writer.writerow(['group_id', 'reference_id', 'event_id', 'event_class', 'included_transcript_id', 'excluded_transcript_id', 'annotated', genotype + '_ref_1', genotype + '_ref_2', genotype + '_ref_3', genotype + '_var_1', genotype + '_var_2', genotype + '_var_3', 'ref_psi', 'var_psi', 'delta_psi', 'chromosome', 'genomic_exon_loc1', 'genomic_exon_loc2', 'strand', 'plasmid_splice_site', 'pval'])
    PSIs = {}
    for event_type in event_types:
        b = open(base_dir + genotype + '_' + event_type + '.psivec', 'r')
        reader = csv.reader(b, delimiter = '\t')
        next(reader)
        for line in reader: 
            PSIs[line[0]] = [float(PSI)*100 for PSI in line[1:]]
        b.close()
        #HepG2_ref_A3-HepG2_var_A3_dPSI	HepG2_ref_A3-HepG2_var_A3_p-val
        #10003;A3:10003:44-219:44-239:+	-0.1698131447	0.0
        d = open(base_dir + genotype + '_' + event_type + '.dpsi', 'r')
        reader = csv.reader(d, delimiter = '\t')
        next(reader)
        for line in reader:
            annotated = False
            this_reference_info = reference_info[line[0]]
            included = this_reference_info['alternative_transcripts']
            excluded = this_reference_info['total_transcripts'].replace(this_reference_info['alternative_transcripts'], '')
            if excluded[0] == ',':
                excluded = excluded[1:]
            if excluded[-1] == ',':
                excluded = excluded[0:-1]
            if 'included' in included and 'excluded' in excluded:
                annotated = True
            new_line = [line[0].split(';')[0], var2ref[line[0].split(';')[0]], line[0], event_type, included, excluded, annotated] + PSIs[line[0]]
            ref_psi = average_psi(PSIs[line[0]][0:3])
            var_psi = average_psi(PSIs[line[0]][3:])
            new_line += [ref_psi, var_psi, float(line[1])*100]
            new_line += add_genomic_locations(new_line) + [line[2]]
            writer.writerow(new_line)
        d.close()
    c.close()



import sys, csv

sys.path.append('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/RBP-amp/rbpamp')
sys.path.append('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/RBP-amp/rbpamp/RBPamp')


import RBPamp
from RBPamp.util import load_model, eval_model, motif_peaks
from collections import defaultdict

RBP_sources = defaultdict(set)
seqs = {}
#chromosome	position	reference	alternative	intron_lower	intron_upper	exon_start	exon_end	strand	control	test_sequencefull_sequence	barcode	group_id
#chr9	114250554	CG	GT	114250541	114250711	114250614	114250668	+	SRSF1_HepG2|8877_SE	CAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTG	CTGACTCTCTCTGCCTCCAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTGCAATTGACTACTAGTGTCCCCTATCTCTAGAGGGCCCGTTTA	GTCCCCTATC	1
a = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/take2/SRE_pool_sequences_dups_hg19.tsv', 'r') 
reader = csv.DictReader(a, delimiter = '\t')
for line in reader:
    RBP_sources[line['group_id']].add(line['control'].split('|')[0])
    seqs[line['group_id']] = line['test_sequence']
a.close()

RBNS_samples = set(['HNRNPK', 'HNRNPL', 'KHSRP', 'RBFOX2', 'RBM22', 'SRSF9', 'TIA1']) 

a = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/NextSeq_Sept_2019/transcriptome_assembled/HepG2_suppa_postprocessed.tsv', 'r')
reader = csv.DictReader(a, delimiter = '\t')
#group_id	reference_id	event_id	event_class	included_transcript_id	excluded_transcript_id	annotated	HepG2_ref_1	HepG2_ref_2	HepG2_ref_3	HepG2_var_1	HepG2_var_2	HepG2_var_3	ref_psi	var_psi	delta_psi	pval
#10000	20401	10000;SE:10000:44-239:283-1003:+	SE	10000_included	10000_excluded	True	0.5555338287055143	0.6393707139975796	0.5357441665891978	0.7681085492756581	0.7999710522506875	0.6428372222738682	0.5768829030974305	0.7369722746000713	0.1600893715	0.018481518500000002
ref_to_ids = defaultdict(set)
ref_to_var = defaultdict(set); var_to_events = defaultdict(set)
for line in reader:
    if line['event_class'] == 'SE':
        unified_id = line['event_id'].split(':')
        unified_id[0] = unified_id[0].replace(line['group_id'], line['reference_id'])
        unified_id[1] = line['reference_id']
        ref_to_ids[line['reference_id']].add(':'.join(unified_id))
        ref_to_var[line['reference_id']].add(line['group_id'])
        var_to_events[line['group_id']].add(':'.join(unified_id))
a.close()

#chromosome	position	reference	alternative	intron_lower	intron_upper	exon_start	exon_end	strand	control	test_sequencefull_sequence	barcode	group_id
#chr9	114250554	CG	GT	114250541	114250711	114250614	114250668	+	SRSF1_HepG2|8877_SE	CAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTG	CTGACTCTCTCTGCCTCCAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTGCAATTGACTACTAGTGTCCCCTATCTCTAGAGGGCCCGTTTA	GTCCCCTATC	1
a = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/take2/SRE_pool_sequences_dups_hg19.tsv', 'r') 
reader = csv.DictReader(a, delimiter = '\t')
seqs = {}; SRE2ref = {}
for line in reader:
    seqs[line['group_id']] = line['test_sequence']
    if line['reference'] == '.':
        SRE2ref[line['control']] = line['group_id']
a.close()


b = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/Variant_RBP-amp/RBPamp_scores.tsv', 'w')
writer = csv.writer(b, delimiter = '\t')
writer.writerow(['rMATs_id', 'group_id', 'reference_id', 'unified_event_id', 'motif_id', 'upstream_ref', 'upstream_var', 'exon_ref', 'exon_var', 'downstream_ref', 'downstream_var'])
for SRE in SRE2ref:
    RBP = SRE.split('_')[0]
    if RBP in RBNS_samples: 
        ref_id = SRE2ref[SRE]
        seq = seqs[ref_id]
        model = load_model('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/RBP-amp/motifs/' + RBP + '.txt')
        offset = int((11 + 1)/2)    
        ref_scores = eval_model(model, seq) 
        for variant in ref_to_var[ref_id]:
            var_seq = seqs[variant]
            var_scores = eval_model(model, var_seq) 
            for splicing_event in var_to_events[variant]:
                three_prime_ss = int(splicing_event.split(':')[2].split('-')[1])
                five_prime_ss = int(splicing_event.split(':')[3].split('-')[0])
                output_line = [SRE, variant, ref_id, splicing_event, RBP]
                output_line += [sum(ref_scores[0:three_prime_ss-157-offset]), sum(var_scores[0:three_prime_ss-157-offset])]
                output_line += [sum(ref_scores[three_prime_ss-157-offset:five_prime_ss-157-offset]), sum(var_scores[three_prime_ss-157-offset:five_prime_ss-157-offset])] 
                output_line += [sum(ref_scores[five_prime_ss-157-offset:]), sum(var_scores[five_prime_ss-157-offset:])]
                writer.writerow(output_line) 
b.close()
#
#

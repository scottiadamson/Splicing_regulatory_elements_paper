import csv, glob, Bio, numpy
from collections import defaultdict
from Bio import motifs

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

all_PWM_files = glob.glob('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/mCross/eCLIP_PWM/*.mat')
RBP_to_PWMs = defaultdict(set); all_PWMs = {};pssm_dict = {}
for PWM_file in all_PWM_files:
    try:
        a = open(PWM_file)
        this_motif = list(motifs.parse(a, "TRANSFAC", strict = False))[0]
        PWM_id = this_motif.get('ID')
        all_PWMs[PWM_id] = this_motif
        RBP = '_'.join(PWM_id.split('.')[0:2][::-1])
        RBP_to_PWMs[RBP].add(PWM_id)
        pssm_dict[PWM_id] = this_motif.pssm
        a.close()
    except IndexError:
        continue 
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


b = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/mCross/mCross_scores.tsv', 'w')
writer = csv.writer(b, delimiter = '\t')
writer.writerow(['rMATs_id', 'group_id', 'reference_id', 'unified_event_id', 'motif_id', 'upstream_ref', 'upstream_var', 'exon_ref', 'exon_var', 'downstream_ref', 'downstream_var'])
for SRE in SRE2ref:
    ref_id = SRE2ref[SRE]
    seq = seqs[ref_id]
    for PWM in RBP_to_PWMs[SRE.split('|')[0]]:
        this_motif = all_PWMs[PWM]
        this_pssm = pssm_dict[PWM]
        ref_scores = numpy.clip(this_pssm.calculate(seq), a_min = 0, a_max = 1000000)
        PWM_offset = int((this_motif.length + 1)/2)
        for variant in ref_to_var[ref_id]:
            var_seq = seqs[variant]
            var_scores = numpy.clip(this_pssm.calculate(var_seq), a_min = 0, a_max = 1000000)
            for splicing_event in var_to_events[variant]:
                three_prime_ss = int(splicing_event.split(':')[2].split('-')[1])
                five_prime_ss = int(splicing_event.split(':')[3].split('-')[0])
                output_line = [SRE, variant, ref_id, splicing_event, PWM]
                output_line += [sum(ref_scores[0:three_prime_ss-157-PWM_offset]), sum(var_scores[0:three_prime_ss-157-PWM_offset])]
                output_line += [sum(ref_scores[three_prime_ss-157-PWM_offset:five_prime_ss-157-PWM_offset]), sum(var_scores[three_prime_ss-157-PWM_offset:five_prime_ss-157-PWM_offset])] 
                output_line += [sum(ref_scores[five_prime_ss-157-PWM_offset:]), sum(var_scores[five_prime_ss-157-PWM_offset:])]
                writer.writerow(output_line) 
b.close()

#10000;SE:10000:44-239:283-1003:+




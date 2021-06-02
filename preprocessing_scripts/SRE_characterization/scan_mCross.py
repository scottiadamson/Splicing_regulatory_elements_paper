import csv, glob, Bio, numpy, gzip
from collections import defaultdict
from Bio import motifs

all_PWM_files = glob.glob('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/mCross/eCLIP_PWM/*.mat')
RBP_to_PWMs = defaultdict(set); all_PWMs = {}; PSSMs = {}
for PWM_file in all_PWM_files:
    try:
        a = open(PWM_file)
        this_motif = list(motifs.parse(a, "TRANSFAC", strict = False))[0]
        PWM_id = this_motif.get('ID')
        all_PWMs[PWM_id] = this_motif
        RBP = '_'.join(PWM_id.split('.')[0:2][::-1])
        RBP_to_PWMs[RBP].add(PWM_id)
        PSSMs[PWM_id] = this_motif.pssm
        a.close()
    except IndexError:
        continue 

#pureclip/HepG2_KHSRP_sensitive_robust.fa 
#>chr1:70697255-70697318(+)_sensitive
#AATGTGCCATGTTGTTTTACCTATCTCTCTTTCTCTCTCACTCCCATGCACACATCCTGTGTG
data_type = 'iCount'

b = gzip.open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/mCross/sensitive_robust/mCross_scores_' + data_type + '.tsv.gz', 'wt')
writer = csv.writer(b, delimiter = '\t')
writer.writerow(['sample', 'status', 'motif_source', 'score'])
fasta_files = glob.glob('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/mCross/sensitive_robust/' + data_type + '/*.fa')
#HepG2_AGGF1_sensitive_robust.fa
for file_ in sorted(fasta_files):
    sample = file_.replace('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/mCross/sensitive_robust/' + data_type + '/', '').replace('_sensitive_robust.fa', '')
    rev_sample = '_'.join(sample.split('_')[::-1])
    if len(RBP_to_PWMs[rev_sample]) == 0:
        if 'K562' in sample:
            if len(RBP_to_PWMs[rev_sample.replace('K562', 'HepG2')]) != 0:
                rev_sample = rev_sample.replace('K562', 'HepG2')
        else:
            if len(RBP_to_PWMs[rev_sample.replace('HepG2', 'K562')]) != 0:
                rev_sample = rev_sample.replace('HepG2', 'K562')
    a = open(file_, 'r')
    print(sample)
    for line in a:
        line = line.strip()
        if line[0] == '>':
            status = line.split('_')[1]  #>chr1:70697255-70697318(+)_sensitive
        else:
            seq = line
            aggregate_scores = []
            for PWM in RBP_to_PWMs[rev_sample]:
                this_motif = all_PWMs[PWM]
                this_pssm = PSSMs[PWM]
                scores = numpy.clip(this_pssm.calculate(seq), a_min = 0, a_max = 1000000)
                aggregate_scores.append(sum(scores)/ (len(scores)- this_motif.length))
            if len(aggregate_scores) != 0:
                writer.writerow([sample, status, '_'.join(rev_sample.split('_')[::-1]), numpy.mean(aggregate_scores)])
    a.close()
b.close()


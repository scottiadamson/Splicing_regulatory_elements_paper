#source activate RBPamp_env
import sys, csv
from scipy import stats

sys.path.append('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/RBP-amp/rbpamp')
sys.path.append('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/RBP-amp/rbpamp/RBPamp')

import RBPamp
from RBPamp.util import load_model, eval_model, motif_peaks

RBNS_samples = 'HepG2_FUBP3,HepG2_PCBP1,K562_FUS,K562_PCBP1,HepG2_FUS,HepG2_PCBP2,K562_HNRNPA1,K562_RBFOX2,HepG2_HNRNPA1,HepG2_RBFOX2,K562_HNRNPC,K562_TARDBP,HepG2_HNRNPC,HepG2_RBM22,K562_HNRNPK,K562_TIA1,HepG2_HNRNPK,HepG2_SRSF9,K562_HNRNPL,K562_TRA2A,HepG2_HNRNPL,HepG2_TIA1,HepG2_KHSRP,HepG2_TRA2A,K562_KHSRP'.split(',')

c = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/sensitive_robust/iCount/RBNS_peak_strength_sensitive_robust.tsv', 'w')
#c = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/sensitive_robust/pureclip/RBNS_peak_strength_sensitive_robust.tsv', 'w')
writer = csv.writer(c, delimiter = '\t')
writer.writerow(['sample', 'peak_id', 'status', 'max_score', 'mean_score'])

peak_scores = {}
for sample in RBNS_samples:
    RBP = sample.split('_')[1]
    a = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/sensitive_robust/iCount/'+ sample + '_sensitive_robust.fa', 'r')
    #a = open('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/sensitive_robust/pureclip/'+ sample + '_sensitive_robust.fa', 'r')
    model = load_model('/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/ENCODE_data/RBNS/RBP-amp/motifs/' + RBP + '.txt')
    RBNS_scores = []
    print(sample)
    for line in a:
        if line[0] == '>':
            seq_id = line.replace('\n', '').replace('>', '')
        else:
            seq = line.replace('\n', '')
            RBNS_score = eval_model(model, seq)
            writer.writerow([sample, seq_id.split('_')[0], seq_id.split('_')[1], max(RBNS_score), sum(RBNS_score)/len(RBNS_score)])
    a.close()
c.close()

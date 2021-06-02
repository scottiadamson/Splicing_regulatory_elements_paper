import csv, sys
from collections import defaultdict
#chr10  100018847   100019486   .   .   -   chr10   100018955   100018956   7.24599;    7.24599 -   HepG2_SRSF7

eCLIP_dir='/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP/'
cell_line = sys.argv[1]
a = open(eCLIP_dir + cell_line + '_rMATs_all_binders.bed', 'r')
reader = csv.reader(a, delimiter = '\t')
exons2RBPs = defaultdict(dict)
for line in reader:
    exon_coordinate = line[0] + ':' + str(int(line[1]) + 300) + '_' + str(int(line[2]) - 300)
    RBP = line[12]
    if RBP not in exons2RBPs[exon_coordinate]:
        exons2RBPs[exon_coordinate][RBP] = 1
    else:
        exons2RBPs[exon_coordinate][RBP] += 1
a.close()

cell_line_samples = {'HepG2':"HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 HepG2_RBM22".split(' '), 'K562': "K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2 K562_FXR1".split(' ')}

sensitive_samples = defaultdict(set)
for sample in cell_line_samples[cell_line]:
    #ID  GeneID  geneSymbol  chr strand  exonStart_0base exonEnd upstreamES  upstreamEE  downstreamES    downstreamEE    ID.1    IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen  SkipFormLen PValue  FDR IncLevel1   IncLevel2   IncLevelDifference
    #1604    ENSG00000227232.4   WASH7P  chr1    -   18267   18366   17914   18061   24737   24891   1604    3,12    3,1 4,4 2,1 197 99  0.87218120344   1.0 0.334,0.858 0.501,0.668 0.011
    a = open(eCLIP_dir + sample + '/SE.MATS.JunctionCountOnly.txt', 'r')
    reader = csv.DictReader(a, delimiter = '\t')
    for line in reader:
        if float(line['FDR']) <=0.05 and abs(float(line['IncLevelDifference'])) >= 0.01:
             exon_coordinate = line['chr'] + ':' + line['exonStart_0base'] + '_' + line['exonEnd']
             sensitive_samples[exon_coordinate].add(sample)
    a.close()

#b = open(eCLIP_dir + cell_line + '_rMATs_binding_summary.bed', 'w')
b = open(eCLIP_dir + cell_line + '_rMATs_binding_summary_iCount.bed', 'w')
writer = csv.writer(b, delimiter = '\t')
writer.writerow(['exon', 'n_RBPs', 'RBPs', 'n_binding_sites', 'sensitive_events'])
for exon in sorted(list(exons2RBPs)):
    writer.writerow([exon, len(exons2RBPs[exon]), ','.join(list(exons2RBPs[exon])), ','.join([str(exons2RBPs[exon][RBP]) for RBP in list(exons2RBPs[exon])]), ','.join(sorted(list(sensitive_samples[exon])))])
b.close()


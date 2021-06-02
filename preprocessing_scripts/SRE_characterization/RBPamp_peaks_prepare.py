#module load anaconda/4.4.0
#source activate pybed
import csv, pybedtools, itertools, sys
from pybedtools import BedTool
from collections import defaultdict

#samples="HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2 HepG2_RBM22 K562_FXR1".split(' ')

samples = ['HepG2_FUBP3', 'HepG2_HNRNPK', 'K562_HNRNPK', 'HepG2_PCBP1', 'K562_PCBP1', 'HepG2_PCBP2', 'HepG2_RBM22', 'HepG2_TIA1', 'K562_TIA1', 'HepG2_FUS', 'K562_FUS', 'HepG2_HNRNPL', 'K562_HNRNPL', 'HepG2_SRSF9', 'HepG2_TRA2A', 'K562_TRA2A', 'HepG2_HNRNPC', 'K562_HNRNPC', 'K562_KHSRP', 'HepG2_KHSRP', 'HepG2_RBFOX2', 'K562_RBFOX2', 'K562_TARDBP', 'HepG2_HNRNPA1', 'K562_HNRNPA1', 'K562_KHDRBS1'] 
#TARDBP_enrichment_R.6mers.tsv

window_size = 30

CLIP_dir = '/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP/'
chrom_sizes = '/home/CAM/adamson/Refs/hg19/chrom.sizes'
genome_fasta = '/home/CAM/adamson/Refs/hg19/GRCh37.p13.genome.fa'

data_type = 'iCount'
#data_type ='pureclip'

for sample in samples:
    print(sample)
    if data_type == 'pureclip':
        CLIP_bed = BedTool(CLIP_dir + sample + '/' + sample + '_eCLIP_binding_regions.bed')
    elif data_type == 'iCount':
        CLIP_rep1 = BedTool(CLIP_dir + sample + '/' + sample + '_rep1_clusters.bed')
        CLIP_rep2 = BedTool(CLIP_dir + sample + '/' + sample + '_rep2_clusters.bed')
        CLIP_bed = CLIP_rep1 + CLIP_rep2.merge()
    CLIP_bed = CLIP_bed.filter(lambda b: b.chrom != 'chrM')
    slop_BT = CLIP_bed.slop(l= window_size ,r = window_size, g = chrom_sizes)
    fasta = pybedtools.example_filename(genome_fasta)

    dPSI_threshold = 0.01; fdr_threshold = 0.05
    #b = open(CLIP_dir + sample + '/' + sample + '_deduped_splicing_output_unnormed.tsv', 'r')
    b = open(CLIP_dir + sample + '/SE.MATS.JunctionCountOnly.txt', 'r')
    reader = csv.DictReader(b, delimiter = '\t')
    KD_sensitive = '';seen_KD = set()
    for line in reader:
        bed_string = ' '.join([line['chr'], str(int(line['exonStart_0base']) - 300), str(int(line['exonEnd']) + 300), line['ID'], '.', line['strand']]) + '\n'
        if float(line['FDR']) <= fdr_threshold and bed_string not in seen_KD:
            if abs(float(line['IncLevelDifference'])) >= dPSI_threshold:
                KD_sensitive += bed_string
                seen_KD.add(bed_string)
    b.close()

    KD_sensitive = pybedtools.BedTool(KD_sensitive, from_string = True)
    sensitive_bed = slop_BT.intersect(KD_sensitive, u = True, s = True)
    robust_bed = slop_BT.intersect(KD_sensitive, v = True, s = True)
    bed_dict = {'sensitive': sensitive_bed, 'robust':robust_bed}
    a = open(CLIP_dir + 'RBNS_peaks/' + data_type + '/' + sample + '_sensitive_robust.fa', 'w')
    for bed_type in ['sensitive', 'robust']:
        bed = bed_dict[bed_type]
        CLIP_seq = bed.sequence(fi = fasta, s = True)
        i = 0
        for line in open(CLIP_seq.seqfn):
            line = line.replace('\n','')
            if line[0] == '>':
                seq_id = line
                a.write(line + '_' + bed_type + '\n')
            else:
                a.write(line + '\n')
            i += 1
        if bed_type == 'sensitive':
            print(i)
    a.close()
    pybedtools.cleanup(remove_all=True)


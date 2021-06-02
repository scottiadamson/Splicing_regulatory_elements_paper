#module load anaconda/4.4.0
#source activate pybed
import csv, pybedtools, itertools, sys, os
from pybedtools import BedTool
from collections import defaultdict

samples="HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2 HepG2_RBM22 K562_FXR1".split(' ')

window_size = 30

CLIP_dir = '/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP/'
chrom_sizes = '/home/CAM/adamson/Refs/hg19/chrom.sizes'
genome_fasta = '/home/CAM/adamson/Refs/hg19/GRCh37.p13.genome.fa'
temp_dir = '/scratch/adamson/tmp2/'
data_type = 'iCount'
#data_type ='pureclip'

c = open('/home/CAM/adamson/eCLIP_peaks_formalized/SRE_number_testable.tsv', 'w')
writer = csv.writer(c, delimiter = '\t')
writer.writerow(['dataset', 'direction', 'total_number', 'testable_number'])
for sample in samples:
    RBP = sample.split('_')[1]
    print(sample)
    if data_type == 'pureclip':
        CLIP_bed = BedTool(CLIP_dir + sample + '/' + sample + '_eCLIP_binding_regions.bed')
    elif data_type == 'iCount':
        CLIP_rep1 = BedTool(CLIP_dir + sample + '/' + sample + '_rep1_clusters.bed')
        CLIP_rep2 = BedTool(CLIP_dir + sample + '/' + sample + '_rep2_clusters.bed')
        CLIP_bed = CLIP_rep1 + CLIP_rep2.merge()
    CLIP_bed = CLIP_bed.filter(lambda b: b.chrom != 'chrM')
    CLIP_bed.saveas(temp_dir + 'test.bed')
    
    dPSI_threshold = 0.01; fdr_threshold = 0.05
    #b = open(CLIP_dir + sample + '/' + sample + '_deduped_splicing_output_unnormed.tsv', 'r')
    b = open(CLIP_dir + sample + '/SE.MATS.JunctionCountOnly.txt', 'r')
    reader = csv.DictReader(b, delimiter = '\t')
    KD_up = ''; KD_down = ''; seen_region = set();testable_up = ''; testable_down = '' 
    for line in reader:
        bed_string = ' '.join([line['chr'], str(int(line['exonStart_0base']) - 300), str(int(line['exonEnd']) + 300), line['ID'], '.', line['strand']]) + '\n'
        if line['strand'] == '+':
            if (int(line['exonEnd']) + 20) - (int(line['exonStart_0base']) - 50) <= 171:
                testable_bed_string = ' '.join([line['chr'], str(int(line['exonStart_0base']) -50), str(int(line['exonEnd']) + 20), line['ID'], '.', line['strand']]) + '\n'
            else:
                testable_bed_string = ''
        else:
            if (int(line['exonEnd']) + 50) - (int(line['exonStart_0base']) - 20) <= 171:
                testable_bed_string = ' '.join([line['chr'], str(int(line['exonStart_0base']) - 20), str(int(line['exonEnd']) + 50), line['ID'], '.', line['strand']]) + '\n'
            else:
                testable_bed_string = ''
        if float(line['FDR']) <= fdr_threshold:
            if float(line['IncLevelDifference']) >= dPSI_threshold:
                if bed_string not in seen_region:
                    KD_up += bed_string
                    testable_up += testable_bed_string
                    seen_region.add(bed_string)
            elif float(line['IncLevelDifference']) <= dPSI_threshold:
                if bed_string not in seen_region:
                    KD_down += bed_string
                    testable_down += testable_bed_string
                    seen_region.add(bed_string)
    b.close()
    KD_up = pybedtools.BedTool(KD_up, from_string = True)
    KD_down = pybedtools.BedTool(KD_down, from_string = True)
    testable_up = pybedtools.BedTool(testable_up, from_string = True)
    testable_down = pybedtools.BedTool(testable_down, from_string = True)
    n_up = len(KD_up.intersect(pybedtools.BedTool(temp_dir + 'test.bed'), u = True, s = True))
    testable_up = len(testable_up.intersect(pybedtools.BedTool(temp_dir + 'test.bed'), u = True, s = True))
    n_down = len(KD_down.intersect(pybedtools.BedTool(temp_dir + 'test.bed'), u = True, s = True))
    testable_down = len(testable_down.intersect(pybedtools.BedTool(temp_dir + 'test.bed'), u = True, s = True))
    writer.writerow([sample, 'up', n_up, testable_up])
    writer.writerow([sample, 'down', n_down, testable_down])
    pybedtools.cleanup(remove_all=True)
    os.remove(temp_dir + 'test.bed')
c.close()

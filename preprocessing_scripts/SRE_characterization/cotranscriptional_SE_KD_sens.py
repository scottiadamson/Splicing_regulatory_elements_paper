#module load anaconda/4.4.0
#source activate pybed
import csv, pybedtools, itertools, sys
import numpy as np
from pybedtools import BedTool
from collections import defaultdict

samples="HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2 HepG2_RBM22 K562_FXR1".split(' ')
pybedtools.helpers.set_tempdir('/scratch/adamson/')

#samples = ['HepG2_FUBP3', 'HepG2_HNRNPK', 'K562_HNRNPK', 'HepG2_PCBP1', 'K562_PCBP1', 'HepG2_PCBP2', 'HepG2_RBM22', 'HepG2_TIA1', 'K562_TIA1', 'HepG2_FUS', 'K562_FUS', 'HepG2_HNRNPL', 'K562_HNRNPL', 'HepG2_SRSF9', 'HepG2_TRA2A', 'K562_TRA2A', 'HepG2_HNRNPC', 'K562_HNRNPC', 'K562_KHSRP', 'HepG2_KHSRP', 'HepG2_RBFOX2', 'K562_RBFOX2', 'K562_TARDBP'] 
#TARDBP_enrichment_R.6mers.tsv
window_size = 30

CLIP_dir = '/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP/'
chrom_sizes = '/home/CAM/adamson/Refs/hg19/chrom.sizes'
genome_fasta = '/home/CAM/adamson/Refs/hg19/GRCh37.p13.genome.fa'
intron_dir = '/home/CAM/adamson/eCLIP_peaks_formalized/motif_expression_regression/'

a = open('/home/CAM/adamson/Refs/hg19/gencode.v19.annotation.gtf', 'r')
reader = csv.reader(a, delimiter = '\t')
gene_to_strand = {}
#chr1   HAVANA  gene    11869   14412   .   +   .   gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
for line in reader:
    if line[0][0] != '#':
        if line[2] == 'gene':
            gene_id = line[8].split(';')[0].replace('gene_id ', '').replace('"', '')
            gene_to_strand[gene_id] = line[6]
a.close()

#Chromosome,Upstream exon region start,Upstream exon region end,Intron start,Intron end,Downstream exon region start,Downstream exon region end,strand,Wild type splicing efficiency,Standard Deviation - 3 replicates,,
b = open(intron_dir + 'Saldi_2021_cotranscriptional_SE_table_S1_updated.csv', 'r')
reader = csv.DictReader(b)
intron_SEs = {}
for line in reader:
    coors = '_'.join([line['Chromosome'], line['Intron start'], str(int(line['Intron end']) -1), line['strand']])
    intron_SEs[coors] = line['Wild type splicing efficiency']
b.close()

for sample in samples:
    print(sample)
    #CLIP_bed = BedTool(CLIP_dir + sample + '/' + sample + '_eCLIP_binding_regions.bed')
    CLIP_rep1 = BedTool(CLIP_dir + sample + '/' + sample + '_rep1_clusters.bed')
    CLIP_rep2 = BedTool(CLIP_dir + sample + '/' + sample + '_rep2_clusters.bed')
    CLIP_bed = CLIP_rep1 + CLIP_rep2.merge()
    CLIP_bed = CLIP_bed.filter(lambda b: b.chrom != 'chrM')

    slop_BT = CLIP_bed.slop(l= window_size ,r = window_size, g = chrom_sizes)

    dPSI_threshold = 0.01; fdr_threshold = 0.05
    #b = open(CLIP_dir + sample + '/' + sample + '_deduped_splicing_output_unnormed.tsv', 'r')
    OG_line = {}
    b = open(CLIP_dir + sample + '/SE.MATS.JunctionCountOnly.txt', 'r')
    reader = csv.DictReader(b, delimiter = '\t')
    KD_sensitive = ''; KD_robust = ''
    for line in reader:
        bed_string = ' '.join([line['chr'], str(int(line['exonStart_0base']) - 300), str(int(line['exonEnd']) + 300), line['ID'], '.', line['strand']]) + '\n'
        if float(line['FDR']) <= fdr_threshold:
            if abs(float(line['IncLevelDifference'])) >= dPSI_threshold:
                KD_sensitive += bed_string
            else:
                KD_robust += bed_string
        else:
            KD_robust += bed_string
        OG_line[line['ID']] = line
    b.close()

    c = open(intron_dir + 'KD_sens_cotranscriptional_SE_iCount/' + sample + '_cotranscriptional_SE.tsv', 'w')
    writer = csv.writer(c, delimiter = '\t')
    writer.writerow(['rMATs_id', 'sensitivity', 'intron_1_coors', 'intron_2_coors', 'skipped_intron_coors', 'intron_1_SE', 'intron_2_SE', 'skipped_SE'])
    KD_sensitive = pybedtools.BedTool(KD_sensitive, from_string = True)
    KD_robust = pybedtools.BedTool(KD_robust, from_string = True)
    sensitive_bed = slop_BT.intersect(KD_sensitive, s = True, wa = True, wb = True)
    robust_bed = slop_BT.intersect(KD_robust, s = True, wa = True, wb = True)
    bed_dict = {'sensitive': sensitive_bed, 'robust':robust_bed}
    peaks = {}; peak_status = {}
    for bed_type in ['sensitive', 'robust']:
        CLIP_bed = bed_dict[bed_type]
        cell_line = sample.split('_')[0]
        seen_ids = set()
        for intersection in CLIP_bed:
            #[u'chr14', u'95085513', u'95085582', u'RP11-986E7.7-ENSG00000273259.1', u'9', u'+', u'chr14', u'95085231', u'95086105', u'7647', u'.', u'+']
            intersection = list(intersection)
            rMATs_id = intersection[9]
            if rMATs_id not in seen_ids:
                rMATs_line = OG_line[rMATs_id]
                #ID      GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE    ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen     PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference
                intron_1_coors = '_'.join([rMATs_line[x] for x in ['chr', 'upstreamEE', 'exonStart_0base', 'strand']])
                intron_2_coors = '_'.join([rMATs_line[x] for x in ['chr', 'exonEnd', 'downstreamES', 'strand']])
                skipped_coors = '_'.join([rMATs_line[x] for x in ['chr', 'upstreamEE', 'downstreamES', 'strand']])
                out_line = [rMATs_id, bed_type, intron_1_coors, intron_2_coors, skipped_coors]
                for coors in [intron_1_coors, intron_2_coors, skipped_coors]:
                    if coors in intron_SEs:
                        out_line.append(intron_SEs[coors])
                    else:
                        out_line.append('NA')
                writer.writerow(out_line)
                seen_ids.add(rMATs_id)
    c.close()
pybedtools.cleanup(remove_all=True)


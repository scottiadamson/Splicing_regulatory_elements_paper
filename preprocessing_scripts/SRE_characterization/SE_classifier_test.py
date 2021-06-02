#module load anaconda2/4.4.0
#source activate pybed

import csv, gzip, pybedtools


samples="HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2 HepG2_RBM22 K562_FXR1".split(' ')

fdr_threshold = 0.05; dPSI_threshold = 0.01
base_dir = '/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP/'
pybedtools.helpers.set_tempdir('/scratch/adamson/tmp2')

def append_bed(event_dict, r1, r2, r3):
    events = [r1, r2, r3]
    for i in range(0,3):
        event_dict[i] += events[i]
        i += 1
    return event_dict

c = open(base_dir + 'SRE_peaks_splicing_events.tsv', 'w')

writer = csv.writer(c, delimiter = '\t')
writer.writerow(['cell_line_RBP', 'region', 'up_peak', 'up_no_peak', 'down_peak', 'down_no_peak', 'background_peak', 'background_no_peak',  'peak_caller'])
seen_up = set(); seen_down = set(); seen_background = set()
for sample in samples:
    print(sample)
        #ID  GeneID  geneSymbol  chr strand  exonStart_0base exonEnd upstreamES  upstreamEE  downstreamES    downstreamEE    ID.1    IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen  SkipFormLen PValue  FDR IncLevel1   IncLevel2   IncLevelDifference
    #1725    ENSG00000227232.4   WASH7P  chr1    -   18267   18366   17914   18061   24737   24891   1725    8,5 2,3 4,4 2,1 197 99  0.911949032935  1.0 0.668,0.456 0.501,0.668 -0.022
    #b = open(base_dir + sample + '/' + sample + '_deduped_splicing_output_unnormed.tsv', 'r')
    b = open(base_dir + sample + '/SE.MATS.JunctionCountOnly.txt', 'r')
    reader = csv.DictReader(b, delimiter = '\t')
    event_type = {'up':['','','',''], 'down':['','','',''], 'background':['','','','']}
    for line in reader:
        if line['strand'] == '+':
            R1 = ' '.join([line['chr'], str(int(line['exonStart_0base']) - 300), str(line['exonStart_0base']), line['ID'], '.', line['strand']]) + '\n' 
            R2 = ' '.join([line['chr'], line['exonStart_0base'], line['exonEnd'], line['ID'], '.', line['strand']]) + '\n'
            R3 = ' '.join([line['chr'], line['exonEnd'], str(int(line['exonEnd']) + 300), line['ID'], '.', line['strand']]) + '\n'
        else:
            R1 = ' '.join([line['chr'], line['exonEnd'], str(int(line['exonEnd']) + 300), line['ID'], '.', line['strand']]) + '\n'
            R2 = ' '.join([line['chr'], line['exonStart_0base'], line['exonEnd'], line['ID'], '.', line['strand']]) + '\n'
            R3 = ' '.join([line['chr'], str(int(line['exonStart_0base']) - 300), line['exonStart_0base'], line['ID'], '.', line['strand']]) + '\n'
        if float(line['FDR']) <= fdr_threshold:
            if float(line['IncLevelDifference']) >= dPSI_threshold and R2 not in seen_up:
                event_type['up'] = append_bed(event_type['up'], R1, R2, R3)
                seen_up.add(R2)
            if float(line['IncLevelDifference']) <= -dPSI_threshold and R2 not in seen_down:
                event_type['down'] = append_bed(event_type['down'], R1, R2, R3)
                seen_down.add(R2)
        elif R2 not in seen_background:
            event_type['background'] = append_bed(event_type['background'], R1, R2, R3)
            seen_background.add(R2)
    b.close()

    #now that we have all the raw count data, let's process the peaks using the various peak callers
    #IDR peaks are just a bed file, easy in
    #chr8   95541294    95541362    SRSF1_HepG2_IDR 1000    -   4.83466847610601    6.64500542862659    -1  -1
    IDR_bed = pybedtools.BedTool(base_dir + sample + '/' + sample + '_IDR_peaks.bed.gz')
    #iCount files need to be joined from the two replicates and foramtted
    iCount_beds = dict()
    for rep in ['1', '2']:
        iCount_beds[rep] = pybedtools.BedTool(base_dir + sample + '/' + sample + '_rep' + rep + '_clusters.bed')
    iCount_bed = iCount_beds['1'].cat(iCount_beds['2'], postmerge = False)
    iCount_bed = iCount_bed.sort()
    iCount_bed = iCount_bed.merge(c = '4,5,6', o = 'distinct,mean,distinct', s = True)
    #Pureclip binding region files are just bed files, so they can just be read in
    #chr1    569897  569898  6.8817; 6.8817  +
    pureclip_bed = pybedtools.BedTool(base_dir + sample + '/' + sample + '_eCLIP_binding_regions.bed')
    pureclip_bc1_bed = pybedtools.BedTool(base_dir + sample + '/' + sample + '_eCLIP_binding_regions_bc1.bed')
    beds = [IDR_bed, iCount_bed, pureclip_bed, pureclip_bc1_bed]
    #Now intersect and do Fisher's test and write output
    j = 0;bed_labels = ['IDR', 'iCount', 'pureclip', 'pureclip_bc1']
    for bed in beds: 
        region_output = {}
        for i in range(0,3):
            region_output[i] = {}
            for etype in ['up', 'down', 'background']:
                region_output[i][etype] = {}
                this_bed = pybedtools.BedTool(event_type[etype][i], from_string = True)
                intersect = this_bed.intersect(bed, s = True, F = 0.5, u =True).count()
                region_output[i][etype]['with'] = intersect
                region_output[i][etype]['without'] = this_bed.count() - intersect
        #['cell_line_RBP', 'region', 'up_peak', 'up_no_peak', 'down_peak', 'down_no_peak', 'background_peak', 'background_no_peak', 'enhancer_FisherExactTest', 'silencer_FisherExactTest']
        region_dict = dict(zip([0,1,2], ['upstream_intron', 'exon', 'downstream_intron']))
        for region in region_dict:
            output_line = [sample, region_dict[region]]
            for etype in ['up', 'down', 'background']:
                output_line += [region_output[region][etype]['with'], region_output[region][etype]['without']]
            writer.writerow(output_line + [bed_labels[j]])
        j += 1
    pybedtools.cleanup(remove_all=True)

c.close()


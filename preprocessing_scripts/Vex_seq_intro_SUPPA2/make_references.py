import csv

upstream_1o = 'CCACTGACTCTCTCTGCCTC'#TGCAG PstI
downstream_1o = 'TCTAGAGGGCCCGTTTAAACCCGCT'#starts with XbaI
upstream_2o = 'AGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTCTAGT'
downstream_2o = 'TCTAGAGGGCCCGTTTAAACCCGCT'
E3 = 'CAATTGTTTTCTTTTGTTTAATTCTTGCTTTCTTTTTTTTTCTTCTCCGCAATTTTTACTATTATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATTAAGTAACTTAAAAAAAAACTTTACACAGTCTGCCTAGTACATTACTATTTGGAATATATGTGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTTTAATTGATACATAATCATTATACATATTTATGGGTTAAAGTGTAATGTTTTAATATCGATACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATTTCTGCATATAAATATTTCTGCATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTCTAGT'#XbaI
RNA_upstream = 'NNGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATATGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGGCCCACTGACTCTCTCTGCCTC'#TGCAG PstI
RNA_downstream = 'TCTAGAGGGCCCGTTTAAACCCGCTGATCAGC' + 'NNNNNNNNNN'
#group_id   test_sequence   full_sequence   barcode
#1   CAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTG  CTGACTCTCTCTGCCTCCAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTGCAATTGACTACTAGTGTCCCCTATCTCTAGAGGGCCCGTTTA   GTCCCCTATC
prefix = '/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/'
a=open('/home/CAM/adamson/Vex-seq/large_pool/SRE_design/SRE_pool_sequence_no_dups.tsv', 'r')
reader = csv.reader(a, delimiter ='\t')
next(reader)
seen = set()
for line in reader:
    b = open(prefix + 'fastas/plasmid/1o/' + str(line[0]) + '_' + str(line[3]) + '.fa', 'w')
    b.write('>' + line[0] + '_' + line[3] + '\n')
    b.write(upstream_1o + line[1] + 'CAATTGACTACTAGT' + line[3] + downstream_1o)
    b.close()
    c = open('fastas/plasmid/2o/' + str(line[0]) + '_' + str(line[3]) + '.fa', 'w')
    c.write('>' + line[0] + '_' + line[3] + '\n')
    c.write(upstream_2o + line[3] + downstream_2o)
    c.close()
a.close() 

print('RNA references')
#chromosome      position        reference       alternative     intron_lower    intron_upper    exon_start      exon_end        strand  control test_sequence   full_sequence   barcode group_id
#chr9    114250554       CG      GT      114250541       114250711       114250614       114250668       +       SRSF1_HepG2|8877_SE     CAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTG      CTGACTCTCTCTGCCTCCAGGGGAAGGAGGTGGAGGGCTGGGTCCCCGGATTCACGTTGTTTCTCTTGCTTTCGGGAATGTGTCTCATAGGGATTTATGGGATTCATTGGTCTGGTCGGGGAGCCAGGAATCGTGGGAGAAAAGGTAAGTGGTGTTGAGGGGAAAAGATAAACAATTAGAGCTTGTGCAATTGACTACTAGTGTCCCCTATCTCTAGAGGGCCCGTTTA   GTCCCCTATC      1
e = open('../SRE_pool_sequences_dups.tsv', 'r')
reader = csv.DictReader(e,delimiter = '\t')
seen = set()
for line in reader:
    if line['group_id'] not in seen:
        f = open('gtfs/' + line['group_id'] + '.gtf', 'w')
        writer = csv.writer(f, delimiter = '\t',quotechar='', escapechar= '\\', quoting=csv.QUOTE_NONE)
        gene_seq = RNA_upstream + line['test_sequence'] + E3 + 'NNNNNNNNNN' + RNA_downstream
        gene_length = len(gene_seq)
        #10000   .       exon    1       46      .       +       .       transcript_id "10000_excluded"; gene_id "10000";
        #10000   .       exon    1005    1197    .       +       .       transcript_id "10000_excluded"; gene_id "10000";
        #10000   .       exon    1       46      .       +       .       transcript_id "10000_included"; gene_id "10000";
        #10000   .       exon    241     285     .       +       .       transcript_id "10000_included"; gene_id "10000"
        #10000   .       exon    1005    1197    .       +       .       transcript_id "10000_included"; gene_id "10000";
        writer.writerow([line['group_id'], '.', 'exon', 1, 46, '.', '+', '.', 'transcript_id "' + line['group_id'] + '_excluded"; gene_id "' + line['group_id'] + '";'])
        if line['strand'] == '+':
            offset = int(line['exon_start']) - int(line['intron_lower'])
            exon_length = int(line['exon_end']) - int(line['exon_start'])
        else:
            offset = int(line['intron_upper']) - int(line['exon_end'])
            exon_length = int(line['exon_end']) - int(line['exon_start'])
        E3_location = gene_seq.find('CTCCTGGGCAACGTGCTGGTCTGTGTGC')+1
        writer.writerow([line['group_id'], '.', 'exon', E3_location, gene_length, '.', '+', '.', 'transcript_id "' + line['group_id'] + '_excluded"; gene_id "' + line['group_id'] + '";'])
        writer.writerow([line['group_id'], '.', 'exon', 1, 46, '.', '+', '.', 'transcript_id "' + line['group_id'] + '_included"; gene_id "' + line['group_id'] + '";'])
        writer.writerow([line['group_id'], '.', 'exon', offset+157, offset + exon_length+156, '.', '+', '.', 'transcript_id "' + line['group_id'] + '_included"; gene_id "' + line['group_id'] + '";'])
        writer.writerow([line['group_id'], '.', 'exon', E3_location, gene_length, '.', '+', '.', 'transcript_id "' + line['group_id'] + '_included"; gene_id "' + line['group_id'] + '";'])
        f.close()
    seen.add(line['group_id'])
e.close()


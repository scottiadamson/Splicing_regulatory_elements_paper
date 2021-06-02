import csv, glob

#1014   1014    44_1003 1014_excluded   216,228,144 216,228,144

base_dir = '/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test/'
genotypes = ['HepG2', 'K562']#, 'K562_NT', 'K562_UPF1']
for genotype in genotypes:
    these_files = glob.glob(base_dir + 'count_files/' + genotype + '*.tsv')
    b = open(base_dir + genotype + '_all_counts.tsv', 'w')
    seen = set()
    for file_ in these_files:
        if genotype == 'K562':
            if 'NT' not in file_ and 'UPF1' not in file_:
                a =open(file_, 'r')
                for line in a:
                    transcript_id = '_'.join(line.split('\t')[2:4])
                    if transcript_id not in seen:
                        b.write(line)
                        seen.add(transcript_id)
                a.close()
        else:
            a =open(file_, 'r')
            for line in a:
                transcript_id = '_'.join(line.split('\t')[2:4])
                if transcript_id not in seen:
                    b.write(line)
                    seen.add(transcript_id)
            a.close()
    b.close()




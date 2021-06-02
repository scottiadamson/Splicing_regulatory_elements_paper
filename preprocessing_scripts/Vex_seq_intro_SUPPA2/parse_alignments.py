import pysam, sys, csv
from collections import defaultdict

variant_id = sys.argv[1]
reference_id = sys.argv[2]
genotype = sys.argv[3]

base_dir = '/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test/' 

def parse_gtf(file_name):
    a = open(base_dir + file_name, 'r')
    reader = csv.reader(a, delimiter = '\t')
    annotated_junctions = {}; junctions = [];junction_names = {}
    for line in reader:
        if line[0][0] != '#':
            if line[2] == 'transcript':
                if junctions != []:
                    junctions = '_'.join([str(x) for x in junctions[1:-1]])
                    annotated_junctions[junctions] = transcript_name
                    junction_names[transcript_name] = junctions
                junctions = []
                transcript_name = line[8].split(';')[1].split(' ')[2].replace('"', '').replace(reference_id, variant_id)
            if line[2] == 'exon':
                junctions.append(int(line[3]) - 1) 
                junctions.append(line[4])
    junctions = '_'.join([str(x) for x in junctions[1:-1]])
    annotated_junctions[junctions] = transcript_name
    junction_names[transcript_name] = junctions
    a.close()
    return annotated_junctions, junction_names

try:
    annotated_junctions, junction_names = parse_gtf('new_gtfs/' + variant_id + '_merged.gtf')
except NameError:
     annotated_junctions, junction_names = parse_gtf('references/' + variant_id +'.gtf')

print(annotated_junctions)
def parse_junction(r):
    BAM_CREF_SKIP = 3 #BAM_CREF_SKIP
    match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
    
    base_position = r.pos
    junctions = []
    for op, nt in r.cigartuples:
        if op in match_or_deletion:
            base_position += nt
        elif op == BAM_CREF_SKIP:
            junc_start = base_position
            base_position += nt
            junctions += [junc_start, base_position]
    return junctions

def parse_bamfile(ID, seq_type, sample):
    if seq_type == 'variant':
        bamfile = base_dir + 'majiq_bams/' + ID + '/' + sample + '_' + ID + '_var_renamed.bam'
    else:
        bamfile = base_dir + 'reference_bams/' + ID + '/' + sample + '_' + ID + '_ref_renamed.bam'
    samfile = pysam.AlignmentFile(bamfile, "rb")
    junction_dict = {}
    new_junctions = defaultdict(lambda:0)
    for read in samfile.fetch():
        if read.is_proper_pair and read.is_secondary == False:
            if read.is_read1:
                junction_dict[read.query_name] = parse_junction(read)
            elif read.query_name in junction_dict:
                read2_junctions = parse_junction(read)
                if read2_junctions != junction_dict[read.query_name]:
                    junction_dict[read.query_name] += parse_junction(read)
                junction = '_'.join([str(x) for x in sorted(list(set(junction_dict[read.query_name])))])
                new_junctions[junction] += 1
    samfile.close()
    return new_junctions

junction_dicts = {'reference': {}, 'variant': {}}
all_junctions = set()
for replicate in ['1', '2', '3']:
    sample = genotype + '_' + replicate
    junction_dicts['variant'][sample] = parse_bamfile(variant_id, 'variant', sample)
    junction_dicts['reference'][sample] = parse_bamfile(reference_id, 'reference', sample)
    all_junctions |= set(junction_dicts['variant'][sample])
    all_junctions |= set(junction_dicts['reference'][sample])

a = open(base_dir + 'count_files/' + genotype + '_' + reference_id + '.tsv', 'a')
writer = csv.writer(a, delimiter = '\t')
#variant_id, control_id, junction, annotation, ref1,ref2,ref3, var1,var2,var3
i = 0
for junction in all_junctions:
    if junction in annotated_junctions:
        if 'cluded' in annotated_junctions[junction]:
            annotation = annotated_junctions[junction]
        else:
            annotation = variant_id + '_novel' + str(i)
        IJC1 = ','.join([str(x) for x in [junction_dicts['reference'][genotype + '_1'][junction], junction_dicts['reference'][genotype + '_2'][junction], junction_dicts['reference'][genotype + '_3'][junction]]])
        IJC2 = ','.join([str(x) for x in [junction_dicts['variant'][genotype + '_1'][junction], junction_dicts['variant'][genotype + '_2'][junction], junction_dicts['variant'][genotype + '_3'][junction]]])
        writer.writerow([variant_id, reference_id, junction, annotation, IJC1, IJC2])
        i+=1
a.close()
 

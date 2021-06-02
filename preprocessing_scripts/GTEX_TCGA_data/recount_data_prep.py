#module load anaconda2/4.4.0
#source activate recount_reg_env

import csv, gzip, numpy, subprocess, os
from collections import defaultdict

base_dir = "/home/CAM/adamson/Recount2_files/"
#this is just an intermediate file paresed from the metadata
#sample_ids.tsv
#0   SRP007461   SRR545720
#71110   TCGA    7794C8A0-3E30-4D9A-B0AE-55EFCAAB8047
#the TCGA ID is the gdc_file_id in the metadata 
a = open(base_dir + 'sample_ids.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
sample_name_dict = {}
for line in reader:
    if line[1] in set(['TCGA', 'SRP012682']):
        sample_name_dict[line[2]] = line[0]
a.close()

#read in junctions of all SREs identified from eCLIP data processing
SRE_junction_dict = {'SE': defaultdict(dict), 'upstream' : defaultdict(dict), 'downstream' : defaultdict(dict)}
b = open('/home/CAM/adamson/eCLIP_peaks_formalized/relevant_SREs_hg38.bed', 'r')
reader = csv.reader(b, delimiter = '\t') 
#chr6    33275796        33275883        AGGF1_HepG2|30835_SE    +
#chr6    33272586        33272726        AGGF1_HepG2|30835_upstream      +
#chr6    33275964        33276066        AGGF1_HepG2|30835_downstream    +
all_junctions = set()
for line in reader:
    junction_type = line[3].split('_')[2]
    junction_id = '_'.join(line[3].split('_')[0:2])
    SRE_junction_dict[junction_type][junction_id] = '_'.join(line[0:3] + [line[4]])
    all_junctions.add('_'.join(line[0:3] + [line[4]]))
b.close()
print(len(all_junctions))

def identify_junction_indicies(project_name):
    a = gzip.open(base_dir + project_name + '/' + project_name + '.junction_id_with_transcripts.bed.gz', 'r')
    reader = csv.reader(a, delimiter = '\t')
    #chr1    10112   10243   4|D:NA|A:NA|J:NA    1000    -
    junction_dict = {}
    for line in reader:
        junction_long_name = '_'.join(line[0:3] + [line[5]])
        junction_id = line[3].split('|')[0]
        if junction_long_name in all_junctions:
            junction_dict[junction_id] = junction_long_name
    a.close()
    return junction_dict

print('Parsing junction ids...')
#####################
#only run this code block (lines 51-66 the first time to write "SRE_relevant_junction_ids.tsv") once, this takes a long time
recount_junctions = identify_junction_indicies('TCGA')
GTEx_junctions = identify_junction_indicies('SRP012682')
##the ones in common are the same ... I checked
recount_junctions.update(GTEx_junctions)
GTEx_junctions = None #just saving memory


print(len(all_junctions))
a = open('/home/CAM/adamson/Recount2_files/SRE_relevant_junction_ids.tsv', 'w')
writer = csv.writer(a, delimiter = '\t')
#junction_id, junction_long_name
for junction_id in recount_junctions:
    if recount_junctions[junction_id] in all_junctions:
        writer.writerow([junction_id, recount_junctions[junction_id]])
a.close()
####################

####################
#If you have already generated "SRE_relevant_junction_ids.tsv", just read it in here and comment out lines 51-66
a = open('/home/CAM/adamson/Recount2_files/SRE_relevant_junction_ids.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
recount_junctions = {}
for line in reader:
    recount_junctions[line[0]] = line[1]
a.close()
####################

def get_junction_counts(project_name):
    a = gzip.open(base_dir + project_name +'/' + project_name + '.junction_coverage.tsv.gz', 'r')
    reader = csv.reader(a, delimiter = '\t')
    #junc_id    sample_id   coverage
    #93606940   68490   1
    #25  64269   2
    junction_count_dict = {}
    for line in reader:
        #get in all_junctions
        junction_id = line[0]
        if junction_id in recount_junctions:
            long_junction_name = recount_junctions[junction_id]
            if long_junction_name in all_junctions:
                sample_ids = line[1].split(',')
                junction_counts = dict(zip(sample_ids, line[2].split(',')))
                junction_count_dict[long_junction_name] = junction_counts
    a.close()
    return junction_count_dict

###############################
#Run this code block (lines 100 - 118) if you haven't written the junction counts to a file yet (generated "SRE_relevant_sample_counts.tsv.gz")
print('Parsing junction counts...')
GTEx_junction_counts = get_junction_counts('SRP012682')
TCGA_junction_counts = get_junction_counts('TCGA')
junction_counts = {}
a = gzip.open('/home/CAM/adamson/Recount2_files/SRE_relevant_sample_counts.tsv.gz', 'w')
writer = csv.writer(a)
junction_counts = {}
for junction_id in recount_junctions:
    long_junction_name = recount_junctions[junction_id]
    output_line = [long_junction_name]; samples = []; counts = []
    for count_dict in [GTEx_junction_counts[long_junction_name], TCGA_junction_counts[long_junction_name]]:
        for sample in count_dict:
            samples.append(sample)
            counts.append(count_dict[sample])
    writer.writerow(output_line + [','.join(samples), ','.join(counts)])
    junction_counts[long_junction_name] = dict(zip(samples, [int(x) for x in counts]))
GTEx_junction_counts = None; TCGA_junction_counts = None
a.close()
###############################

###############################
#If you've alredy generated "SRE_relevant_sample_counts.tsv.gz" just run this block (lines 121-127) and comment out the first block (lines 100 - 118) 
a = gzip.open('/home/CAM/adamson/Recount2_files/SRE_relevant_sample_counts.tsv.gz', 'rt')
reader = csv.reader(a, delimiter = '\t')
junction_counts = {}
for line in reader:
    junction_counts[line[0]] = dict(zip(line[1].split(','), [int(x) for x in line[2].split(',')]))
a.close()
###############################

def legit_line(line, project_name):
    if project_name == 'TCGA':
        if len(line) == 859:
            return True
        else:
            return False
    else:
        return True

def parse_meta_data(project_name):
    #I want to filter out any samples with potential GoF (missense) mutations in RBPs
    if project_name == 'TCGA':
        TCGA_keepers = set()
        a = open(base_dir + 'TCGA_no_RBP_mut_case_ids.txt', 'r')
        #ff38628e-8abe-44d6-8d3b-7cfe46bd2e35 gdc_cases.case_id
        #unfortuntely TCGAbiolinks returns the case_id, while recount is using the file_id as an identifier
        for line in a:
            TCGA_keepers.add(line.replace('\n', ''))
        a.close()
    field_codes = {'tissue': {'SRP012682': 'smts', 'TCGA' : 'gdc_cases.project.project_id'}, 'batch' : {'SRP012682': 'smgebtch', 'TCGA': 'cgc_case_batch_number'}, 'sample_id': {'SRP012682': 'run', 'TCGA': 'gdc_file_id'}}
    meta_data_dict = {}
    a = open(base_dir + project_name + '/' + project_name + '.tsv', 'r')
    reader = csv.reader(a, delimiter= '\t')
    header = next(reader)
    for line in reader:
        #for whatever reason, some of the lines have random line endings in the middle of them... I'm just going to filter them out
        if legit_line(line, project_name):
            line = dict(zip(header, line))
            sample_number = sample_name_dict[line[field_codes['sample_id'][project_name]].upper()]
            retain_sample = True
            if project_name == 'TCGA':
                if line['gdc_cases.case_id'] not in TCGA_keepers:
                    retain_sample = False
            if retain_sample:
                meta_data_dict[sample_number] = {}
                for field in field_codes:
                    value = line[field_codes[field][project_name]]
                    if value == None:
                        value = '0'
                    meta_data_dict[sample_number][field] = value.upper()
                for field in ['avg_read_length', 'paired_end']:
                    meta_data_dict[sample_number][field] = line[field]
    a.close()
    return meta_data_dict

GTEx_metadata = parse_meta_data('SRP012682')
metadata = parse_meta_data('TCGA')
metadata.update(GTEx_metadata)
GTEx_metadata = None

#calculate PSI for each sample
PSI_dict = {}; weight_dict = {}
for SRE_id in SRE_junction_dict['SE']:
    excluded_junction = SRE_junction_dict['SE'][SRE_id]
    upstream_inc_junction = SRE_junction_dict['upstream'][SRE_id]
    downstream_inc_junction = SRE_junction_dict['downstream'][SRE_id]
#    #check that these juncitons are all observed in the external datasets
    legit_junctions = 0
    for junction_type in [excluded_junction, upstream_inc_junction, downstream_inc_junction]:
        if junction_type in junction_counts:
            legit_junctions += 1
    included_count_dict = {}; excluded_count_dict = {}
    if legit_junctions == 3:
        all_samples = set(junction_counts[excluded_junction]) | set(junction_counts[upstream_inc_junction]) | set(junction_counts[downstream_inc_junction])
        for sample in all_samples:
            for dict_ in [upstream_inc_junction, downstream_inc_junction, excluded_junction]:
                if sample not in junction_counts[dict_]:
                    junction_counts[dict_][sample] = 0
            excluded_count_dict[sample] = junction_counts[excluded_junction][sample] 
            included_count_dict[sample] = (junction_counts[upstream_inc_junction][sample] + junction_counts[downstream_inc_junction][sample])
        exon_start_coor = int(upstream_inc_junction.split('_')[2])
        exon_end_coor = int(downstream_inc_junction.split('_')[1])
        exon_length = exon_end_coor - exon_start_coor
        PSI_values = []; weights = [];samples = []
        for sample in all_samples:
            try:
                # I was going to do this with the rMATs formula, but the read length is unavailable from TCGA samples
                # I'm just going to use the inc junction counts / (inc junction counts + skipped junction counts*2)
                weight = included_count_dict[sample] + excluded_count_dict[sample]
                if weight >=10 and sample in metadata:
                    PSI = included_count_dict[sample] / (included_count_dict[sample] + (excluded_count_dict[sample]*2))
                    #don't trust anything with less than a weight of 10
                    PSI_values.append(round(PSI, 3))
                    weights.append(weight)
                    samples.append(sample)
            except KeyError:
                continue #don't add the sample to anything
        PSI_dict[SRE_id] = dict(zip(samples, PSI_values))
        weight_dict[SRE_id] = dict(zip(samples, weights))

#read in RBPs for input variables
RBP_dir = '/home/CAM/adamson/Multi-CRISPRs/'
RBP_files = ['Hentze_2018_RBPs.csv', 'Gerstberger_2014_RBPs.csv']
RBPs = set() #formatted as versionless ENSEMBL gene Ids
for file_ in RBP_files:
    a = open(RBP_dir + file_, 'r')
    reader = csv.DictReader(a)
    for line in reader:
        if file_ == 'Hentze_2018_RBPs.csv':
            for RBP in line['ID'].split(';'):
                RBPs.add(RBP)
        else:
            RBPs.add(line['gene id'])
RBPs.remove("")
print(len(RBPs))

#get gene info
a = open(base_dir + 'gene_info.tsv', 'r')
reader = csv.DictReader(a, delimiter = '\t')
#gene_id    bp_length
#ENSG00000000003.14  4535
gene_lengths = {}
for line in reader:
    gene_lengths[line['gene_id'].split('.')[0]] = int(line['bp_length'])
a.close()


#read in genes
def import_gene_expression(project_name):
    a = gzip.open(base_dir + project_name + '/counts_gene.tsv.gz', 'r')
    reader = csv.DictReader(a, delimiter = '\t')
    #for TCGA these are all CAPS case_ids, and for GTEx they are SRR ids
    #I downloaded the v2 files, so the gene_id is the last column with a versioned gene_id
    samples = set();sample_totals = {}
    for sample in metadata:
        samples.add(metadata[sample]['sample_id'].upper())
        sample_totals[sample] = 0
    gene_counts = defaultdict(dict)
    for line in reader:
        for sample in line:
            gene_id = line['gene_id'].split('.')[0]
            if sample in samples:
                sample_number = sample_name_dict[sample]
                scaled_count = float(line[sample]) / gene_lengths[gene_id]
                sample_totals[sample_number] += scaled_count
                if gene_id in RBPs:
                    gene_counts[gene_id][sample_number] = scaled_count
    a.close()
    for sample in samples:
        sample_number = sample_name_dict[sample]
        if sample_totals[sample_number] != 0:
            #TPM transform here
            for gene_id in gene_counts:
                gene_counts[gene_id][sample_number] = gene_counts[gene_id][sample_number]*1000000/sample_totals[sample_number]
    return gene_counts

##################################
#only run this code block (lines 277-306) if you haven't generated "Recount_RBP_TPMs.tsv.gz" yet, otherwise comment it out and read in the file (line 309-318)
GTEx_TPMs = import_gene_expression('SRP012682')
TCGA_TPMs = import_gene_expression('TCGA')
TPMs = {}
for RBP in RBPs:
    TPMs[RBP] = {}
    for sample in GTEx_TPMs[RBP]:
        TPMs[RBP][sample] = GTEx_TPMs[RBP][sample]
    for sample in TCGA_TPMs[RBP]:
        TPMs[RBP][sample] = TCGA_TPMs[RBP][sample]
GTEx_TPMs = None; TCGA_TPMs = None

measured_samples = set()
for RBP in TPMs:
    for sample in TPMs[RBP]:
        measured_samples.add(sample)
make sure to get all samples

a = gzip.open(base_dir + 'Recount_RBP_TPMs.tsv.gz', 'w')
writer = csv.writer(a, delimiter = '\t')
writer.writerow(['gene_id'] + sorted(list(measured_samples)))
for RBP in sorted(RBPs):
    output_line = [RBP]
    for sample in sorted(list(measured_samples)):
        if sample in TPMs[RBP]:
            output_line.append(TPMs[RBP][sample])
        else:
            output_line.append('NA') #there are lots of NAs here becasue of old gene identifiers mostly
    writer.writerow(output_line)
a.close()
##################################
#If you've already generated "Recount_RBP_TPMs.tsv.gz", comment out lines 277-306
a = gzip.open(base_dir + 'Recount_RBP_TPMs.tsv.gz', 'rt')
reader = csv.DictReader(a, delimiter = '\t')
for line in reader:
    #there are lots of NAs here becasue of old gene identifiers mostly
    #they show up as lines with all NAs; so just don't read them in
    if line['50099'] != 'NA':
        TPMs[line['gene_id']] = line
a.close()
print(len(PSI_dict))


for SRE_id in list(PSI_dict)[0:118]:
    if SRE_id not in already_done:
        b = open(base_dir + 'regression_input_tables/' + SRE_id.replace('|', '_') + '.tsv', 'wt')
        writer = csv.writer(b, delimiter = '\t')
        gene_vals = {}
        for gene in TPMs:
            gene_mean = numpy.mean([float(TPMs[gene][sample]) for sample in PSI_dict[SRE_id]])
            gene_std = numpy.std([float(TPMs[gene][sample]) for sample in PSI_dict[SRE_id]])
            if gene_mean != 0.0 and gene_std != 0.0:
                gene_vals[gene] = {'mean': gene_mean, 'std': gene_std}
        writer.writerow(['sample_id', 'PSI', 'weight', 'batch', 'tissue', 'cancer'] + sorted(list(gene_vals)))
        for sample in sorted(PSI_dict[SRE_id]):
            if 'TCGA' in metadata[sample]['tissue']:
                cancer = 1
            else:
                cancer = 0
            output_line = [sample, PSI_dict[SRE_id][sample], weight_dict[SRE_id][sample], metadata[sample]['batch'], metadata[sample]['tissue'], cancer]
            for gene in sorted(gene_vals):
                output_line.append((float(TPMs[gene][sample]) - gene_mean) / gene_std)
            writer.writerow(output_line)
        b.close() 
        #If you've already done the elastic net regression, comment out line 342, otherwise, run that first (with 343 commented out), and then repeat with selected variables and 343 running
        #regression_command = 'Rscript ' + base_dir + 'logistic_regression.R --sample ' + SRE_id.replace('|', '_')
        regression_command = 'Rscript ' + base_dir + 'logistic_regression_selected_vars.R --sample ' + SRE_id.replace('|', '_')
        print(regression_command.split(' '))
        subprocess.call(regression_command.split(' '))
        os.remove(base_dir + 'regression_input_tables/' + SRE_id.replace('|', '_') + '.tsv')


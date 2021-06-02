#'''GET the results of a search from an ENCODE server'''
#source activate web_parsers

import csv, subprocess, os

base_dir ='/home/CAM/adamson/eCLIP_peaks_formalized/'
cell_lines = ['K562', 'HepG2']
#hardcode U2AFs for controls
desired_samples = {'K562':set(['U2AF1', 'U2AF2']), 'HepG2':set(['U2AF1', 'U2AF2'])}
#RBP,"Essential Genes","Splicing regulation",Spliceosome,"Novel RBP",Other,"eCLIP HepG2","eCLIP K562","RNAseq HepG2","RNAseq K562",RBNS
for cell_line in cell_lines:
    a = open(base_dir + cell_line + '_eCLIP_RNA.csv', 'r')
    reader = csv.DictReader(a)
    for line in reader:
        if line['Spliceosome'] == '0':
            if line['Splicing regulation'] == '1' or line['Novel RBP'] == '1':
                 desired_samples[cell_line].add(line['RBP'])
    a.close()

def download_files(accession, out_file, extension):
    download_command = ['wget', '-q', '-O', base_dir + out_file, 'https://www.encodeproject.org/files/' + accession + '/@@download/' + accession + extension]
    #print(' '.join(download_command))
    #subprocess.call(download_command)

#for splicing data
#RBP,Cell_line,RNA-Seq_exp,Control_exp,RBP knockdown fastq files,RBP knockdown bam files,RBP knockdown RSEM,RBP knockdown CUFFDIFF,RBP knockdown DESeq,RBP knockdown MISO,RBP knockdown rMATS,RBP knockdown DESeq after batch Correction,RBP knockdown rMATS after batch Correction,Non-target control datasets,Control fastq files,Control bam files,Control RSEM
#STAU1,HepG2,ENCSR124KCF,ENCSR067GHD,"ENCFF756OXT,ENCFF036XQB,ENCFF084DHZ,ENCFF691HLO","ENCFF501MSW,ENCFF365QAH","ENCFF834FOW,ENCFF408GMM",ENCFF523JIV,ENCFF436JDD,ENCFF378LQU,ENCFF727RBG,ENCFF577UMY,ENCFF907FID,non-target,"ENCFF773HLT,ENCFF768JSC,ENCFF914FKR,ENCFF838TBG","ENCFF774FAK,ENCFF103PVI","ENCFF990KGR,ENCFF992TTT"
desired_samples['HepG2'].add('RBM22')
desired_samples['K562'].add('FXR1')
#From Van Nostrand et al., 2020 (https://doi.org/10.1038) Supplementary Data 2
b = open(base_dir + 'ENCODE_paper_SD_Table_2_KD-RNA-seq.csv', 'r')
reader = csv.DictReader(b)
d = open(base_dir + 'ENCODE_file_accessions.tsv', 'w')
writer = csv.writer(d, delimiter = '\t')
writer.writerow(['RBP', 'cell_line', 'file_type', 'sample_type', 'file_accession'])
for line in reader:
    if line['RBP'] in desired_samples[line['Cell_line']]:
        download_files(line['RBP knockdown rMATS after batch Correction'], 'Splicing_batch_normed/' + line['Cell_line'] + '_' + line['RBP'] + '.tsv.tar.gz', '.tar.gz')
        download_files(line['RBP knockdown rMATS'], 'Splicing_unnormed/' + line['Cell_line'] + '_' + line['RBP'] + '.tsv.tar.gz', '.tar.gz')
        writer.writerow([line['RBP'], line['Cell_line'], 'rMATS', 'all_RNA-seq_samples', line['RBP knockdown rMATS']])
b.close()

#clip data 
#RBP,Cell_line,eCLIP_exp,Control_exp,Experiment Fastq Files,Experiment Bam Files,Peak Files,Reproducible Peak Files,Input ID,Input Fastq Files,Input Bam Files
#HNRNPC,HepG2,ENCSR550DVK,ENCSR497ANA,"ENCFF035QUK,ENCFF584MZO,ENCFF040BGS,ENCFF065BIC","ENCFF869PNM,ENCFF736ZQQ","ENCFF401SPQ,ENCFF738FJS",ENCFF105AKX,HNRNPC eCLIP mock input,"ENCFF412YDY,ENCFF357SNM",ENCFF823ZZZ
c = open(base_dir + 'ENCODE_paper_SD_Table_2_eCLIP.csv', 'r')
reader = csv.DictReader(c)
for line in reader:
    if line['RBP'] in desired_samples[line['Cell_line']]:
        #os.mkdir(base_dir + 'eCLIP/' + line['Cell_line'] + '_' + line['RBP'])
        download_files(line['Input Bam Files'], 'eCLIP/' + line['Cell_line'] + '_' + line['RBP'] + '_input.bam', '.bam') 
        download_files(line['Experiment Bam Files'].split(',')[0], 'eCLIP/' + line['Cell_line'] + '_' + line['RBP'] + '_rep1.bam', '.bam') 
        download_files(line['Experiment Bam Files'].split(',')[1], 'eCLIP/' + line['Cell_line'] + '_' + line['RBP'] + '_rep2.bam', '.bam') 
        download_files(line['Peak Files'].split(',')[0], 'eCLIP/' +  line['Cell_line'] + '_' + line['RBP'] + '_rep1_peaks.bed.gz', '.bed.gz')
        download_files(line['Peak Files'].split(',')[1], 'eCLIP/' +  line['Cell_line'] + '_' + line['RBP'] + '_rep2_peaks.bed.gz', '.bed.gz')
        download_files(line['Reproducible Peak Files'], 'eCLIP/' +  line['Cell_line'] + '_' + line['RBP'] + '_IDR_peaks.bed.gz', '.bed.gz')
        writer.writerow([line['RBP'], line['Cell_line'], 'eCLIP_bam', 'input', line['Input Bam Files']])
        writer.writerow([line['RBP'], line['Cell_line'], 'eCLIP_bam', 'rep1', line['Experiment Bam Files'].split(',')[0]])
        writer.writerow([line['RBP'], line['Cell_line'], 'eCLIP_bam', 'rep2', line['Experiment Bam Files'].split(',')[1]])
c.close();d.close()




#!/bin/bash
#SBATCH --job-name=umi_tools
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --mem=10G
#SBATCH -o umi_tools.o%j
#SBATCH -e umi_tools.e%j
#SBATCH --qos=general
#SBATCH --mail-type=END

module load anaconda/4.4.0
source activate umi_root
export PYTHONPATH=""

base_dir=/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019
#samples="HepG2_1 HepG2_2 HepG2_3 K562_1 K562_2 K562_3"
#for sample in $samples; do echo $sample; sbatch --export=sample=$sample /home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/mark_umi_duplicates.sh; done 
echo $sample
umi_tools extract -I $base_dir/$sample"_R2.fastq.gz" --bc-pattern=NNNNNNNNNN --read2-in=$base_dir/$sample"_R1.fastq.gz" --stdout=$base_dir/$sample"_R2_processed.fastq.gz" --read2-out=$base_dir/$sample"_R1_processed.fastq.gz" --extract-method=string
#rm $base_dir/$sample"_R2.fastq.gz"
#rm $base_dir/$sample"_R1.fastq.gz"


#!/bin/bash
#SBATCH --job-name=VS_plas_map
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --mem=5G
#SBATCH -o VS_plasmid_map.o%j
#SBATCH -e VS_plasmid_map.e%j
#SBATCH --qos=general

module load bwa/0.7.17
module load samtools
module load bcftools/1.10.2
base_dir="/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307"
refs=$(<$base_dir/barcode_samples.txt)

echo "1o"
for ref in $refs; do
    bwa index $base_dir/fastas/plasmid/1o/$ref".fa"
    bwa mem -t 4 -I 240,24,265,225 $base_dir/fastas/plasmid/1o/$ref".fa" $base_dir/sorted_fastqs/1o_SRE/$ref"_R1.fastq" $base_dir/sorted_fastqs/1o_SRE/$ref"_R2.fastq" > $base_dir/alignments/1o_SRE/$ref".sam"
    samtools view -b -h -@ 4 -F 0x800 -f 3 $base_dir/alignments/1o_SRE/$ref".sam" > $base_dir/alignments/1o_SRE/$ref".bam"
    samtools sort -@ 4 $base_dir/alignments/1o_SRE/$ref".bam" > $base_dir/alignments/1o_SRE/$ref"_sorted.bam"
    samtools mpileup -u -v --output-tags DP,AD --max-depth 10000 -f $base_dir/fastas/plasmid/1o/$ref".fa" $base_dir/alignments/1o_SRE/$ref"_sorted.bam" | bcftools call - -mv -O v -o $base_dir/vcfs/1o_SRE/$ref".vcf"
    echo $ref $(samtools view -@ 4 -F 0x800 -c -f 67 $base_dir/alignments/1o_SRE/$ref"_sorted.bam") >> $base_dir/1o_uniquely_assigned.txt
    cat $base_dir/vcfs/1o_SRE/$ref".vcf" | grep -v "#" >> $base_dir/vcfs/1o_merged.vcf2
    rm $base_dir/alignments/1o_SRE/$ref".sam"
    rm $base_dir/fastas/plasmid/1o/$ref".fa."*
    rm $base_dir/alignments/1o_SRE/$ref".bam"
    rm $base_dir/alignments/1o_SRE/$ref"_sorted.bam"
    rm $base_dir/vcfs/1o_SRE/$ref".vcf"
done
echo "2o"

for ref in $refs; do
    bwa index $base_dir/fastas/plasmid/2o/$ref".fa"
    bwa mem -t 4 -I 298,29,315,275 $base_dir/fastas/plasmid/2o/$ref".fa" $base_dir/sorted_fastqs/2o_SRE/$ref"_R1.fastq" $base_dir/sorted_fastqs/2o_SRE/$ref"_R2.fastq" > $base_dir/alignments/2o_SRE/$ref".sam"
    samtools view -b -h -@ 4 $base_dir/alignments/2o_SRE/$ref".sam" > $base_dir/alignments/2o_SRE/$ref".bam"
    samtools sort -@ 4 $base_dir/alignments/2o_SRE/$ref".bam" > $base_dir/alignments/2o_SRE/$ref"_sorted.bam"
    echo $ref $(samtools view -@ 4 -c -f 67 -F 0x800 $base_dir/alignments/2o_SRE/$ref"_sorted.bam") >> $base_dir/2o_uniquely_assigned.txt
    cat $base_dir/vcfs/2o_SRE/$ref".vcf" | grep -v "#" >> $base_dir/vcfs/2o_merged.vcf
    rm $base_dir/vcfs/2o_SRE/$ref".vcf"
    rm $base_dir/alignments/2o_SRE/$ref".sam"
    rm $base_dir/fastas/plasmid/2o/$ref".fa."*
    rm $base_dir/alignments/2o_SRE/$ref".bam"
    rm $base_dir/alignments/2o_SRE/$ref"_sorted.bam"
done

python /home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_190307/plasmid_summary.py

for ref in $refs; do
    rm $base_dir/sorted_fastqs/1o_SRE/$ref*.fastq
    rm $base_dir/sorted_fastqs/2o_SRE/$ref*.fastq
done


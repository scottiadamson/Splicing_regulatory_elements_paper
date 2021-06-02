#!/bin/bash
#SBATCH --job-name=pureclip_call
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=35G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH -o pureclip_call.o%j
#SBATCH -e pureclip_call.e%j

module load gcc
module load samtools
module load anaconda2/4.4.0

base_dir="/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP"
ref_fasta="/home/CAM/adamson/Refs/hg19/GRCh37.p13.genome.fa"
samples="HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2 K562_FXR1"
export TMPDIR=/scratch/adamson
export WINEXTRACT=/home/CAM/adamson/.conda/envs/pureclip_env/bin/winextract

source activate umi_root
source deactivate
unset PYTHONPATH

for sample in $samples; do
    echo $sample
    R1=$sample"_rep1"
    R2=$sample"_rep2"
    control=$sample"_input"
    files="$R1 $R2 $control"
    source activate umi_root

    for file in $files; do
        samtools view -h $base_dir/$sample/$file".bam" | python $base_dir/fix_UMI_tags.py pureclip - > $base_dir/$sample/$file"_pureclip.sam"
        samtools view -hb -@ 7 $base_dir/$sample/$file"_pureclip.sam" > $base_dir/$sample/$file"_pureclip.bam"
        samtools index $base_dir/$sample/$file"_pureclip.bam"
        umi_tools dedup -I $base_dir/$sample/$file"_pureclip.bam" --paired -S $base_dir/$sample/$file"_pureclip_dedup.bam"
        samtools view -hb -@ 7 -f 130 $base_dir/$sample/$file"_pureclip_dedup.bam" > $base_dir/$sample/$file"_pureclip_dedup_R2.bam"
        samtools index $base_dir/$sample/$file"_pureclip_dedup_R2.bam"
        rm $base_dir/$sample/$file"_pureclip.sam"
        rm $base_dir/$sample/$file"_pureclip.bam"
        rm $base_dir/$sample/$file"_pureclip_dedup.bam"
    done

    source deactivate
    source activate pureclip_env

    compute_CLmotif_scores.sh $ref_fasta $base_dir/$sample/$sample"_rep1_pureclip_dedup_R2.bam" $base_dir/motifs.xml $base_dir/motifs.txt $base_dir/$sample/fimo_clmotif_occurences_rep1.bed
    compute_CLmotif_scores.sh $ref_fasta $base_dir/$sample/$sample"_rep2_pureclip_dedup_R2.bam" $base_dir/motifs.xml $base_dir/motifs.txt $base_dir/$sample/fimo_clmotif_occurences_rep2.bed
    cat $base_dir/$sample/fimo_clmotif_occurences_rep1.bed $base_dir/$sample/fimo_clmotif_occurences_rep2.bed | sort -k1,1 -k2,2n - > $base_dir/$sample/fimo_merged_sorted.bed
    bedtools merge -s -c 4,5,6 -o max,max,distinct -i $base_dir/$sample/fimo_merged_sorted.bed > $base_dir/$sample/fimo_merged_max.bed
    pureclip -i $base_dir/$sample/$sample"_rep1_pureclip_dedup_R2.bam" -bai $base_dir/$sample/$sample"_rep1_pureclip_dedup_R2.bam.bai" -i $base_dir/$sample/$sample"_rep2_pureclip_dedup_R2.bam" -bai $base_dir/$sample/$sample"_rep2_pureclip_dedup_R2.bam.bai" -g $ref_fasta -o $base_dir/$sample/$sample"_eCLIP_CL_peaks.bed" -or $base_dir/$sample/$sample"_eCLIP_binding_regions.bed" -nt 8 -iv 'chr1;chr2;chr3' -ibam $base_dir/$sample/$sample"_input_pureclip_dedup_R2.bam" -ibai $base_dir/$sample/$sample"_input_pureclip_dedup_R2.bam.bai" -fis $base_dir/$sample/fimo_merged_max.bed
    pureclip -i $base_dir/$sample/$sample"_rep1_pureclip_dedup_R2.bam" -bai $base_dir/$sample/$sample"_rep1_pureclip_dedup_R2.bam.bai" -i $base_dir/$sample/$sample"_rep2_pureclip_dedup_R2.bam" -bai $base_dir/$sample/$sample"_rep2_pureclip_dedup_R2.bam.bai" -g $ref_fasta -o $base_dir/$sample/$sample"_eCLIP_CL_peaks_bc1.bed" -or $base_dir/$sample/$sample"_eCLIP_binding_regions_bc1.bed" -nt 8 -iv 'chr1;chr2;chr3' -ibam $base_dir/$sample/$sample"_input_pureclip_dedup_R2.bam" -ibai $base_dir/$sample/$sample"_input_pureclip_dedup_R2.bam.bai" -fis $base_dir/$sample/fimo_merged_max.bed -bc 1
    rm $base_dir/$sample/$sample"_input_pureclip_dedup_R2.bam.bai"
    rm $base_dir/$sample/$sample"_input_pureclip_dedup_R2.bam"
    rm $base_dir/$sample/$sample"_rep2_pureclip_dedup_R2.bam.bai"
    rm $base_dir/$sample/$sample"_rep2_pureclip_dedup_R2.bam"
    rm $base_dir/$sample/$sample"_rep1_pureclip_dedup_R2.bam.bai"
    rm $base_dir/$sample/$sample"_rep1_pureclip_dedup_R2.bam"
    rm $base_dir/$sample/fimo_clmotif_occurences_rep1.bed
    rm $base_dir/$sample/fimo_clmotif_occurences_rep2.bed
    rm $base_dir/$sample/fimo_merged_sorted.bed
    rm $base_dir/$sample/fimo_merged_max.bed
done

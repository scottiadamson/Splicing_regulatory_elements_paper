#!/bin/bash
#SBATCH --job-name=iCount
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --mem=35G
#SBATCH -o iCount.o%j
#SBATCH -e iCount.e%j
#SBATCH --qos=general
#SBATCH --mail-type=END

base_dir="/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP"
temp_dir="/scratch/adamson"
samples="HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2"
module load anaconda/4.4.0
source activate iCount_env

unset PYTHONPATH

export TMP_ROOT="/local/tmp/"

gtf="/home/CAM/adamson/Refs/hg19/gencode.v19.annotation.gtf"
for sample in $samples; do 
    echo $sample
    for rep in {1..2}; do
        echo $sample"_rep"$rep
        export ICOUNT_TMP_ROOT="/local/tmp/iCount/"$sample
        samtools view -h $base_dir/$sample/$sample"_rep"$rep".bam" | python $base_dir/fix_UMI_tags.py iCount - > $base_dir/$sample/$sample"_rep"$rep"_iCount.sam"
        samtools view -hb -@ 7 $base_dir/$sample/$sample"_rep"$rep"_iCount.sam" > $base_dir/$sample/$sample"_rep"$rep"_iCount.bam"
        iCount xlsites $base_dir/$sample/$sample"_rep"$rep"_iCount.bam" $base_dir/$sample/$sample"_rep"$rep"_unique.bed" $base_dir/$sample/$sample"_rep"$rep"_multi.bed" $base_dir/$sample/$sample"_rep"$rep"_skipped.bam" --group_by start --quant cDNA
        rm $base_dir/$sample/$sample"_rep"$rep"_iCount.bam"
        rm $base_dir/$sample/$sample"_rep"$rep"_iCount.bam.bai"
        rm $base_dir/$sample/$sample"_rep"$rep"_iCount.sam"
        rm -r $ICOUNT_TMP_ROOT
        iCount peaks --scores $base_dir/$sample/$sample"_rep"$rep"_scores.tsv" $gtf $base_dir/$sample/$sample"_rep"$rep"_unique.bed" $base_dir/$sample/$sample"_rep"$rep"_peaks.bed"
        rm -r $ICOUNT_TMP_ROOT
        iCount clusters $base_dir/$sample/$sample"_rep"$rep"_unique.bed" $base_dir/$sample/$sample"_rep"$rep"_peaks.bed" $base_dir/$sample/$sample"_rep"$rep"_clusters.bed"
        rm -r $ICOUNT_TMP_ROOT
    done
done


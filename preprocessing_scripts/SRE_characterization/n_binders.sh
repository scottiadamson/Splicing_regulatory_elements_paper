module load bedtools

chrom_sizes="/home/CAM/adamson/Refs/hg19/chrom.sizes"

#HepG2_samples="HepG2_AGGF1 HepG2_BCCIP HepG2_DDX55 HepG2_DDX59 HepG2_FAM120A HepG2_FUBP3 HepG2_FUS HepG2_GRWD1 HepG2_GTF2F1 HepG2_HLTF HepG2_HNRNPA1 HepG2_HNRNPC HepG2_HNRNPK HepG2_HNRNPL HepG2_HNRNPM HepG2_HNRNPU HepG2_KHSRP HepG2_MATR3 HepG2_NCBP2 HepG2_NKRF HepG2_PCBP1 HepG2_PCBP2 HepG2_PPIG HepG2_PTBP1 HepG2_QKI HepG2_RBFOX2 HepG2_RBM15 HepG2_SDAD1 HepG2_SLTM HepG2_SND1 HepG2_SRSF1 HepG2_SRSF7 HepG2_SRSF9 HepG2_SUB1 HepG2_SUGP2 HepG2_TBRG4 HepG2_TIA1 HepG2_TIAL1 HepG2_TRA2A HepG2_U2AF1 HepG2_U2AF2 HepG2_UCHL5 HepG2_XRCC6 HepG2_RBM22"

K562_samples="K562_AATF K562_AGGF1 K562_AKAP8L K562_DDX55 K562_EWSR1 K562_FAM120A K562_FMR1 K562_FUS K562_GRWD1 K562_GTF2F1 K562_HLTF K562_HNRNPA1 K562_HNRNPC K562_HNRNPK K562_HNRNPL K562_HNRNPM K562_HNRNPU K562_KHDRBS1 K562_KHSRP K562_MATR3 K562_METAP2 K562_NCBP2 K562_NONO K562_PCBP1 K562_PHF6 K562_PPIL4 K562_PTBP1 K562_QKI K562_RBFOX2 K562_RBM15 K562_SAFB2 K562_SAFB K562_SDAD1 K562_SLTM K562_SND1 K562_SRSF1 K562_SRSF7 K562_TARDBP K562_TBRG4 K562_TIA1 K562_TRA2A K562_U2AF1 K562_U2AF2 K562_UCHL5 K562_UTP3 K562_WRN K562_XRCC6 K562_ZC3H8 K562_ZNF800 K562_ZRANB2 K562_FXR1"

eCLIP_dir="/home/CAM/adamson/eCLIP_peaks_formalized/eCLIP"
temp_dir="/scratch/adamson"
touch $temp_dir/K562_exons.bed
for sample in $K562_samples; do
    echo $sample
    #chr1   568027  568028  11.0262;    11.0262 +
    touch $temp_dir/K562_eCLIP_peaks.bed
    #awk -v sample="$sample" 'BEGIN{OFS="\t";} {print $0,sample}' $eCLIP_dir/$sample/$sample"_eCLIP_binding_regions.bed" |grep -v "chrM" >> $temp_dir/K562_eCLIP_peaks.bed
    for rep in {1..2}; do 
        awk -v sample="$sample" 'BEGIN{OFS="\t";} {print $0,sample}' $eCLIP_dir/$sample/$sample"_rep"$rep"_clusters.bed" |grep -v "chrM"| sort -k1,1 -k2,2n -  > $temp_dir/K562_$sample"_eCLIP_peaks.bed"
    done
    bedtools merge -i $temp_dir/K562_$sample"_eCLIP_peaks.bed" -s -c 4,5,6,7 -o collapse,collapse,collapse,collapse >> $temp_dir/K562_eCLIP_peaks.bed
    rm $temp_dir/K562_$sample"_eCLIP_peaks.bed"
    awk 'BEGIN{OFS="\t";} {print $4,$6,$7,".",".",$5}' $eCLIP_dir/$sample/$sample"_deduped_splicing_output_unnormed.tsv" |grep -v "exonStart_0base" >> $temp_dir/K562_exons.bed 
done

sort -k 1,1 -k2,2n $temp_dir/K562_exons.bed |sort | uniq > $temp_dir/K562_unique_exons.bed
rm $temp_dir/K562_exons.bed
bedtools slop -i $temp_dir/K562_unique_exons.bed -g $chrom_sizes -b 300  |  sort -k 1,1 -k2,2n - > $temp_dir/K562_unique_exons_slopped.bed
sort -k 1,1 -k2,2n $temp_dir/K562_eCLIP_peaks.bed > $temp_dir/K562_eCLIP_peaks_sorted.bed
bedtools intersect -sorted -s -a $temp_dir/K562_unique_exons_slopped.bed -b $temp_dir/K562_eCLIP_peaks_sorted.bed -wa -wb > $eCLIP_dir/K562_rMATs_all_binders.bed 
python $eCLIP_dir/parse_rMATs_binders.py K562

#rm $temp_dir/*.bed

#!/bin/bash
#SBATCH --job-name=VS_hisat_remap
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=5G
#SBATCH -o hisat_remap.o%j
#SBATCH -e hisat_remap.e%j

module load samtools
base_dir="/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test"
hisat_dir="/isg/shared/apps/hisat2/2.1.0"

base_dir="/home/CAM/adamson/Vex-seq/large_pool/SRE_design/NextSeq_Sept_2019/majiq_test"
samples="HepG2_1 HepG2_2 HepG2_3 K562_1 K562_2 K562_3"
genotypes="HepG2 K562"
#note that this was run piecemeal (a few references at a time) otherwise it will take forever
refs="3333 4882 9246 9404 10882 14437 16724 8144 21820 10417 33292 28895 18615 24441 1192 19313 27200 366 6483 22919 30615 24017 21572 22798 4101 29322 32534 24999 18403 13565 17181 3758 2437 26446 3743 22174 6410 16855 1684 29237 26378 9496 4433 1987 852 8840 22334 29422 7218 8546 23272 2368 10485 15458 20375 28873 19926 26972 13913 24338 19490 22912 6688 32632 12410 16395 29291 18162 15288 30306 21478 26332 6394 2601 27687 25989 28384 14607 13380 19390 27561 32481 28701 7483 10408 20978 31842 18881 9613 29899 16311 17689 17749 11358 22223 25281 14406 25564 13118 31465 20401 20279 20679 20142 1836 2596 22355 17269 4705 8542 15937 16296 28937 16104 5371 19971 7099 30428 15036 2845 11192 3719 3483 25905 13279 30103 27957 258 28888 12577 20209 8225 1420 14030 28769 10725 18111 26856 12959 32761 32983 4081 12693 23634 27042 32022 30038 25900 26501 10874 3147 30977 19880 30968 24976 29517 23536 19526 323 28313 8595 1014 22017 4463 19118 11511 13476 25534 23841 20011 20600 18205 2924 32472 8154 15650 14395 24134 7511 3666 17483 23637 15038 4163 5172 24576 26693 29608 30044 27389 398 10145 26703 3930 7277 3543 24819 23196 11893 16703 6823 9442 2702 549 3506 54 21421 30120 2180 15578 26374 21227 4333 22508 14104 23687 19067 4233 23090 15843 14858 33231 2081 3387 27690 1900 26553 31177 21519 23465 1853 21597 28061 5515 14003 29471 14872 5675 7835 33222 31630 28323 16519 2352 28748 12562 17564 14612 4123 15780 14571 4286 685 15217 1602 28045 27147 22790 23133 29746 19533 9540 3817 23600 28413 25020 25757 21760 7023 2354 16470 16599 11733 29158 23539 7280 21552 479 28327 21970 18755 5970 19562 5495 5345 25971 20936 17397 23617 28046 13393 26076 16137 11151 21586 33136"
for ref in $refs; do
    echo $ref
    mkdir /scratch/adamson/$ref
    TMPDIR=/scratch/adamson/$ref
    readarray variants < $base_dir/control_variant_arrays/$ref".txt"
    module load hisat2/2.1.0
    for variant in $variants; do
        echo $variant
        python $base_dir/make_reference_files.py $variant
        python $hisat_dir/hisat2_extract_splice_sites.py $base_dir/references/$variant".gtf" > $base_dir/references/$variant"_ss.tsv"
        python $hisat_dir/hisat2_extract_exons.py $base_dir/references/$variant".gtf" > $base_dir/references/$variant"_exons.tsv"
        hisat2-build --ss $base_dir/references/$variant"_ss.tsv" --exon $base_dir/references/$variant"_exons.tsv" $base_dir/references/$variant".fa" $base_dir/references/$variant
        mkdir $base_dir/majiq_bams/$variant
    done
    for sample in $samples; do
        echo $sample
        python $base_dir/sort_fastqs2.py $sample $ref
        for variant in $variants; do
            echo $variant
            mkdir $base_dir/references/$variant
            hisat2 --no-softclip --dta-cufflinks --fr -x $base_dir/references/$variant -1 $base_dir/temp_fastqs/$sample/$variant"_R1.fastq.gz" -2 $base_dir/temp_fastqs/$sample/$variant"_R2.fastq.gz" | samtools view -hb - > $base_dir/majiq_bams/$variant/$sample"_var.bam" 
            samtools sort $base_dir/majiq_bams/$variant/$sample"_var.bam" > $base_dir/majiq_bams/$variant/$sample"_var_sorted.bam"
            samtools index $base_dir/majiq_bams/$variant/$sample"_var_sorted.bam"
            rm $base_dir/majiq_bams/$variant/$sample"_var.bam"
        done
    done
    module rm hisat2
    module load anaconda/4.4.0
    source activate umi_root
    export PYTHONPATH=""
    for sample in $samples; do
        for variant in $variants; do
            umi_tools dedup -I $base_dir/majiq_bams/$variant/$sample"_var_sorted.bam" --paired -S $base_dir/majiq_bams/$variant/$sample"_var_dedup.bam" --log2stderr --chimeric-pair discard --unpaired-reads discard 
            rm $base_dir/majiq_bams/$variant/$sample"_var_sorted.bam"
            mv $base_dir/majiq_bams/$variant/$sample"_var_dedup.bam" $base_dir/majiq_bams/$variant/$sample"_var.bam"
            samtools index $base_dir/majiq_bams/$variant/$sample"_var.bam"
        done
    done
    source deactivate
    module load stringtie/2.0.3
    for variant in $variants; do
        echo $variant
        #formely here python $base_dir/change_reference.py $variant ...
        samtools view -H $base_dir/majiq_bams/$ref/"HepG2_1_var.bam" > $base_dir/references/$variant/header.sam
        rm $base_dir/references/$variant*".ht2" 
        for sample in $samples; do
            samtools view $base_dir/majiq_bams/$variant/$sample"_var.bam" | awk -v ref=$ref 'BEGIN {OFS="\t"} {$3=ref} {print $0}' - | cat $base_dir/references/$variant/header.sam - | samtools view -hb - > $base_dir/majiq_bams/$variant/$sample"_"$variant"_var_renamed.bam" 
            samtools index $base_dir/majiq_bams/$variant/$sample"_"$variant"_var_renamed.bam"
            stringtie $base_dir/majiq_bams/$variant/$sample"_"$variant"_var_renamed.bam" -t --fr -G $base_dir/references/$ref".gtf" -o $base_dir/references/$sample"_"$variant"_var_renamed.gtf"
        done
        stringtie --merge -G $base_dir/references/$ref".gtf" -o $base_dir/references/$variant"_merged.gtf" $base_dir/references/*"_"$variant"_var_renamed.gtf" $base_dir/references/*"_"$ref"_var_renamed.gtf"
        cp $base_dir/references/$variant"_merged.gtf" $base_dir/new_gtfs/$variant"_merged.gtf"
    done
    module rm stringtie
    mkdir $base_dir/reference_bams/$ref
    source activate pysam_env
    for sample in $samples; do
        cp $base_dir/majiq_bams/$ref/$sample"_"$ref"_var_renamed.bam" $base_dir/reference_bams/$ref/$sample"_"$ref"_ref_renamed.bam"
        cp $base_dir/majiq_bams/$ref/$sample"_"$ref"_var_renamed.bam.bai" $base_dir/reference_bams/$ref/$sample"_"$ref"_ref_renamed.bam.bai"
    done
    for genotype in $genotypes; do
        touch $base_dir/count_files/$genotype"_"$ref".tsv"
    done
    for variant in $variants; do
        for genotype in $genotypes; do
            python $base_dir/parse_alignments.py $variant $ref $genotype 
        done
        rm -r $base_dir/temp_fastqs/*/$variant"_R"*.fastq.gz
        rm -r $base_dir/majiq_bams/$variant
        rm $base_dir/references/$variant"_"*
        rm $base_dir/references/$variant"."*
        rm -r $base_dir/references/$variant
        rm $base_dir/references/*"_"$variant"_"*.gtf
        rm $base_dir/new_gtfs/$variant"_merged.gtf"
    done
    rm $base_dir/reference_bams/$ref/*_var_renamed.bam
    rm $base_dir/reference_bams/$ref/*_var_renamed.bam.bai
    rm -r $base_dir/reference_bams/$ref
    rm -r /scratch/adamson/$ref
    source deactivate
    module rm anaconda
done
##

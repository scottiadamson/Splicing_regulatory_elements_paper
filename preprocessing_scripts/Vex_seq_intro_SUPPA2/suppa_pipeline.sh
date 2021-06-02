#source activate suppa2_env
suppa_PATH="/Users/sadamson/miniconda3/envs/suppa2_env/bin"
base_dir="/Users/sadamson/Desktop/scripts_and_such/Vex-seq/RBP_pool/NextSeq_Sept_2019/transcriptome_assembled"

python $base_dir/harmonize_and_classify.py
python $suppa_PATH/suppa.py generateEvents -i $base_dir/master_novel_merged.gtf -o $base_dir/master_novel_merged -f ioe -e SE SS MX FL -t 1

genotypes="HepG2 K562" # K562_NT K562_UPF1"
##genotypes="HepG2"
variant_types="ref var"
event_types="A3 A5 SE"
##event_types="A3"
for genotype in $genotypes; do
	for event_type in $event_types; do
		for variant_type in $variant_types; do
			python $suppa_PATH/suppa.py psiPerEvent -e $base_dir/$genotype"_all_counts_suppa_"$variant_type".tsv" -i $base_dir/master_novel_merged_$event_type"_strict.ioe" -o $base_dir/$genotype"_"$variant_type"_"$event_type
		done
		python $suppa_PATH/suppa.py diffSplice --method empirical --input $base_dir/master_novel_merged_$event_type"_strict.ioe" --psi $base_dir/$genotype"_ref_"$event_type.psi $base_dir/$genotype"_var_"$event_type.psi --tpm $base_dir/$genotype"_all_counts_suppa_ref.tsv" $base_dir/$genotype"_all_counts_suppa_var.tsv" --output $base_dir/$genotype"_"$event_type
	done
done	

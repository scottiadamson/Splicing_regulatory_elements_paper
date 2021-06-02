These are the input files, scripts and output files that were used to create data analyzed in the jupyter notebooks.  
Note that the purpose of this is to publically show the analysis that corresponds to the paper _NOT_ to readily use in other contexts without modification first.  
Note that this also generates very large files.  

|input file(s) | script | output | notes |
|---|---|---|---|
| | download_recount2_data.R | gene expression and junction count tables for GTEx and TCGA | |
| | ID_TCGA_mutations.R | TCGA_no_RBP_mut_case_ids.txt | |
| TCGA and GTEx metadata files | make_sample_ids.py | sample_ids.tsv | | 
| gene expression and junction count tables for GTEx and TCGA, gene_info.tsv, TCGA_no_RBP_mut_case_ids.txt, Hentze_2018_RBPs.csv, Gerstberger_2014_RBPs.csv | recount_data_prep.py | regression_outputs_trad/, regression_outputs/, regression_input_tables/ | also runs logistic_regression.R and logistic_regression_selected_vars.R from this script | 
| regression_outputs_trad/ | aggregate_n_p.py | n_p_reg_summary.tsv | | 

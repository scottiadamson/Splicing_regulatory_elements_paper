These are the input files, scripts and output files that were used to create data analyzed in the jupyter notebooks.  
Note that the purpose of this is to publically show the analysis that corresponds to the paper _NOT_ to readily use in other contexts without modification first.

|input file(s) | script | output |
|---|---|---|
| ENCODE_paper_SD_Table_2_KD-RNA-seq.csv, ENCODE_paper_SD_Table_2_KD-RNA-seq.csv | fetch_files.py | eCLIP bam files, ENCODE_file_accessions.tsv | 
| | download_rMATS.sh | rMATs SE files |
| | fetch_files.py | eCLIP bam files |
| eCLIP bam files | iCount.sh | eCLIP postprocessed files |
| eCLIP bam files | pureclip_call.sh | eCLIP postprocessed files | 
| rMATs SE files, eCLIP postprocessed files | SE_classifier_test.py | SRE_peaks_splicing_events.tsv |
| rMATs SE files, eCLIP postprocessed files | RBPamp_peaks_prepare.py | sensitve and robust peak fasta files|
| sensitve and robust peak fasta files | RBNS_scan.py | RBNS_peak_strength_sensitive_robust.tsv |
| sensitve and robust peak fasta files | scan_mCross.py | mCross_scores_iCount.tsv.gz |
| rMATs SE files, eCLIP postprocessed files | parse_rMATs_binders.py, n_binders.sh | HepG2_rMATs_binding_summary_iCount.bed, K562_rMATs_binding_summary_iCount.bed |
| rMATs SE files, eCLIP postprocessed files, Saldi_2021_cotranscriptional_SE_table_S1_updated.csv | cotranscriptional_SE_KD_sens.py | KD_sens_cotranscriptional_SE_iCount/|
| rMATs SE files, eCLIP postprocessed files | SRE_number_testable.py | SRE_number_testable.tsv |

These are the input files, scripts and output files that were used to create data analyzed in the jupyter notebooks.  
Note that the purpose of this is to publically show the analysis that corresponds to the paper _NOT_ to readily use in other contexts without modification first.

|input file(s) | script | output | notes |
|---|---|---|---|
| SRE_pool_sequence_no_dups.tsv | make_references.py | plasmid and RNA-seq fasta files, gtf files | outputs can be concatenated into fasta and gtf files with all references "all.fa", "all.gtf" |
| plasmid fastq files, reference fasta files | map_plasmids.sh, plasmid_summary.py | plasmid_statistics2.tsv | 1o_plasmid_VAFs.tsv and 2o_plasmid_quantifications.tsv are just selected columns from this table|
| raw RNA-seq fastq files | mark_umi_duplicates.sh | marked RNA-seq fastq files | |
| marked RNA-seq fastq files | sort_fastqs_initial.py | initially sorted RNA-seq fastqs | |
| initially sorted RNA-seq fastq files, plasmid_statistics2.tsv, fasta and gtf references | hisat_map.sh, make_reference_files.py, sort_fastqs2.py, parse_alignments.py | tsv count files | |
| tsv count files | summarize_count_files.py | sample summarized count files | |
| all.gtf, sample summarized count files | suppa_pipeline.sh, harmonize_and_classify.py | raw suppa output | |  
| raw suppa output | post_process_suppa.py | K562_suppa_postprocessed.tsv, HepG2_suppa_postprocessed.tsv | | 

library(recount)

download_study('SRP012682', type = 'counts-gene', outdir = 'GTEX')
download_study('TCGA', type = 'counts-gene')

download_study('SRP012682', type = 'counts-exon', outdir = 'GTEX')
download_study('TCGA', type = 'counts-exon')

download_study('SRP012682', type = 'counts-jx', outdir = 'GTEX')
download_study('TCGA', type = 'counts-jx')


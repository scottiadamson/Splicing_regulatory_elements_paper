library(dplyr)
library(TCGAbiolinks)
library(tidyverse)
projects <- getGDCprojects()$project_id

TCGA_projects <- projects[grepl("TCGA", projects)]
split_projects <- strsplit(TCGA_projects, '-', fixed = TRUE)
TCGA_projects <- unlist(lapply(split_projects, function(l) l[[2]]))

RBP_dir <- '/home/CAM/adamson/Multi-CRISPRs/'
Henzte_RBPs <- read.csv(paste0(RBP_dir,'Hentze_2018_RBPs.csv'))
Gerstberger_RBPs <- read.csv(paste0(RBP_dir, 'Gerstberger_2014_RBPs.csv'))
Henzte_RBPs <- unlist(strsplit(as.character(Henzte_RBPs$ID), ';', fixed = TRUE))
Gerstberger_RBPs <- Gerstberger_RBPs$gene.id

all_RBPs <- union(Gerstberger_RBPs, Henzte_RBPs)

input_samples <- 0
retained_samples <- c()
for (cancer in TCGA_projects){
    maf <- GDCquery_Maf(cancer, pipelines = "mutect")
    maf <- maf[,c('Hugo_Symbol', 'Variant_Classification', 'case_id', 'Transcript_ID', 'Gene')]
    input_samples = input_samples + length(unique(maf$Tumor_Sample_Barcode))
    maf <- filter(maf, !(Gene %in% all_RBPs), Variant_Classification %in% c('Nonsense_Mutation', 'Missense_Mutation'))
    retained_samples <- c(retained_samples, unique(maf$case_id))
}

print(paste0('input tumor sample barcodes: ', as.character(input_samples)))
print(paste0('retained tumor sample barcodes: ', as.character(length(retained_samples))))
base_dir <- '/home/CAM/adamson/Recount2_files'
write(retained_samples, paste0(base_dir, '/TCGA_no_RBP_mut_case_ids.txt'))


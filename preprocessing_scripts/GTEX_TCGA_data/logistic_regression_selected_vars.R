library(tidyverse)
library(glmnet)
library(car)
library(boot)
library(R.utils)

sample <- cmdArg("sample")
print(sample)

base_dir <- '/home/CAM/adamson/Recount2_files/'

input_df <- read_tsv(paste0(base_dir, 'regression_input_tables/', sample, '.tsv'))

input_df$PSI <- car::logit(input_df$PSI, adjust = 0.01)
input_df$batch <- as.factor(input_df$batch)
input_df$tissue <- as.factor(input_df$tissue)
input_df$batch <- fct_explicit_na(input_df$batch, na_level = "Unknown")
input_df$tissue <- fct_explicit_na(input_df$tissue, na_level = "Unknown")
input_df$cancer <- as.factor(input_df$cancer)

selected_variables <- read_tsv(paste0(base_dir, 'regression_outputs/', sample, '_coefficients.tsv'), col_names = c('variable', 'coefficient'))
filter_out_variables <- pull(selected_variables %>% filter(coefficient == 0, grepl('ENSG', variable)) %>% dplyr::select(variable))

#AKAP8L_K562_22255_coefficients.tsv
filter_out <- c('sample_id', 'weight', filter_out_variables)
gene_df <- select(input_df, -filter_out)
logistic_model <- glm(formula = paste0('PSI ~ ', paste0(colnames(gene_df)[2:ncol(gene_df)], collapse = ' + ')), weights = input_df$weight, data = gene_df)
res <- as.data.frame(coef(summary(logistic_model)))
names(res) <- c('coefficient', 'std_error', 't_val', 'p_val')
res$variable_id <- rownames(res)
res <- as_tibble(res)
write_tsv(res, paste0(base_dir, 'regression_outputs_trad/', sample, '_coefficients.tsv'))
predictions <- as.data.frame(cbind(input_df$sample_id, inv.logit(input_df$PSI), as.vector(inv.logit(logistic_model$fitted.values))))
colnames(predictions) <- c('sample_id', 'observed_PSI', 'predicted_PSI')

write_tsv(predictions, paste0(base_dir, 'regression_outputs_trad/', sample, '_predictions.tsv'))


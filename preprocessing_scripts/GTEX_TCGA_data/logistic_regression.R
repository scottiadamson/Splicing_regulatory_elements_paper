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

filter_out <- c('sample_id', 'PSI', 'weight')
gene_df <- select(input_df, -filter_out)

reg_mat <- model.matrix( ~ ., gene_df)
logistic_model <- cv.glmnet(reg_mat, input_df$PSI, weights = input_df$weight, alpha = 0.5, maxit = 1000000)

coef.apprx = as.data.frame(as.matrix(coef(logistic_model, s = "lambda.min", exact = FALSE))) 

write.table(coef.apprx, paste0(base_dir, 'regression_outputs/', sample, '_coefficients.tsv'), sep = '\t', row.names = TRUE, quote=FALSE, col.names =FALSE)
predictions <- predict(logistic_model, newx = reg_mat)
predictions <- cbind(input_df$sample_id, inv.logit(input_df$PSI), as.vector(inv.logit(predictions)))
colnames(predictions) <- c('sample_id', 'observed_PSI', 'predicted_PSI')
write.table(predictions, paste0(base_dir, 'regression_outputs/', sample, '_predictions.tsv'), , sep = '\t', row.names = FALSE, quote = FALSE)

jpeg(paste0(base_dir, 'regression_plots/', sample, '_MSE.jpg'))
plot(logistic_model)


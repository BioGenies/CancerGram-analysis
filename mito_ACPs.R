library(drake)
library(biogram)
library(dplyr)
library(ranger)

source("functions/get_mers.R")
source("functions/mc_model_functions.R")
source("functions/train_model_peptides.R")

loadd(peptide_model_mc_anticp)
loadd(imp_ngrams_mc_anticp)
loadd(mer_model_mc_anticp)
peptides <- read_fasta("data/mito_ACPs.fasta")
mers <- mer_df_from_list(peptides)
ngrams <- count_imp_ngrams(mers, imp_ngrams_mc_anticp)

mer_preds <- cbind(mers, predict(mer_model_mc_anticp, as.matrix(ngrams[, imp_ngrams_mc_anticp]))[["predictions"]])
stats <- calculate_statistics_mc(mer_preds,  c("acp", "amp", "neg"))

peptide_preds <- cbind(stats[, c("source_peptide")],
                       predict(peptide_model_mc_anticp,
                               select(stats, -source_peptide))[["predictions"]])
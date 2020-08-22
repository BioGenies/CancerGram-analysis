library(drake)
library(biogram)
library(dplyr)
library(ranger)

source("functions/get_mers.R")
source("functions/mc_model_functions.R")
source("functions/train_model_peptides.R")

loadd(c(peptide_model_mc_anticp, imp_ngrams_mc_anticp, mer_model_mc_anticp, 
        imp_ngrams_mc, mer_model_mc, peptide_model_mc), 
      path = "benchmark_cache")

peptides <- read_fasta("data/mito_ACPs.fasta")
mers <- mer_df_from_list(peptides)

predict_mito_ACPs <- function(mers, imp_ngrams, mer_model, peptide_model) {
  ngrams <- count_imp_ngrams(mers, imp_ngrams)
  mer_preds <- cbind(mers, predict(mer_model, as.matrix(ngrams[, imp_ngrams]))[["predictions"]])
  stats <- calculate_statistics_mc(mer_preds,  c("acp", "amp", "neg"))
  peptide_preds <- cbind(stats[, c("source_peptide")],
                         predict(peptide_model, stats)[["predictions"]])
  return(list("mer_predictions" = mer_preds, "peptide_predictions" = peptide_preds))
}

our_datasets_res <- predict_mito_ACPs(mers, imp_ngrams_mc, mer_model_mc, peptide_model_mc)
anticp_datasets_res <- predict_mito_ACPs(mers, imp_ngrams_mc_anticp, mer_model_mc_anticp, peptide_model_mc_anticp)


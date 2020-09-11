library(dplyr)
library(drake)
library(biogram)
library(magrittr)
library(seqR)
library(ranger)
library(ggplot2)
library(tidyr)
library(cvTools)
library(measures)
library(pbapply)
library(seqinr)
library(xtable)

source("functions/mc_model_functions.R")
source("functions/do_cv.R")
source("functions/train_model_peptides.R")
source("functions/benchmark_functions.R")
source("functions/publication_functions.R")
source("functions/get_mers.R")


if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/Projekty/BioNgramProjects/CancerGram/publication_files/"
}
  
# This drake plan requires targets generated in AntiCP_benchmark.R script by running 
# benchmark_first_models drake plan. 

loadd(c(mer_df_mc, ngrams_mc, mers_mc_anticp, ngrams_mc_anticp,
        benchmark_peptide_preds_mc_anticp, pos_train_main, pos_test_main, 
        neg_train_main, neg_test_main, neg_train_alt, neg_test_alt,
        imp_ngrams_mc_anticp, mer_model_mc_anticp, peptide_model_mc_anticp), 
      path = "benchmark_cache")


publication_results <- drake_plan(
  cv_mer_mc_anticp = do_cv_mc(mers_mc_anticp, ngrams_mc_anticp, 0.05),
  cv_peptide_mc_anticp = do_cv_peptides_mc(cv_mer_mc_anticp),
  cv_mer_mc_anticp_0.001 = do_cv_mc(mers_mc_anticp, ngrams_mc_anticp, 0.001),
  cv_peptide_mc_anticp_0.001 = do_cv_peptides_mc(cv_mer_mc_anticp_0.001),
  cv_mer_performance_measures_0.05 = calc_cv_performance(cv_mer_mc_anticp),
  cv_mer_performance_measures_0.001 = calc_cv_performance(cv_mer_mc_anticp_0.001),
  cv_peptide_perfromance_measures_0.05 = calc_cv_performance(cv_peptide_mc_anticp),
  cv_peptide_perfromance_measures_0.001 = calc_cv_performance(cv_peptide_mc_anticp_0.001),  
  cv_table = get_cv_pred_table(cv_mer_performance_measures_0.05, cv_peptide_perfromance_measures_0.05),
  datasets_table = get_datasets_table(pos_train_main, pos_test_main, neg_train_main, neg_test_main, neg_train_alt, neg_test_alt),
  mito_ACP_preds = predict_mito_ACPs("data/mito_ACPs.fasta", imp_ngrams_mc_anticp, 
                                     mer_model_mc_anticp, peptide_model_mc_anticp),
  mito_ACP_table = get_mito_ACP_pred_table(mito_ACP_preds)
)

make(publication_results, seed = 2938, cache = new_cache("publication_cache"))
  
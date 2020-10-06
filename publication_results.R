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
  
# This drake plan requires targets generated in analysis.R script by running 
# analysis_CancerGram drake plan. 

loadd(c(mer_df_mc, ngrams, imp_ngrams, imp_ngrams_dat, mer_model, 
        peptide_model, benchmark_mer_preds, benchmark_peptide_preds,
        pos_train_main, pos_test_main, neg_train_main, 
        neg_test_main, neg_train_alt, neg_test_alt))


publication_results <- drake_plan(
  cv_mer = do_cv_mc(mers_df_mc, ngrams, 0.0001),
  cv_peptide = do_cv_peptides_mc(cv_mer),
  cv_mer_performance_measures = calc_cv_performance(cv_mer),
  cv_peptide_performance_measures = calc_cv_performance(cv_peptide),
  cv_table = get_cv_pred_table(cv_mer_performance_measures, cv_peptide_performance_measures),
  datasets_table = get_datasets_table(pos_train_main, pos_test_main, neg_train_main, neg_test_main, neg_train_alt, neg_test_alt),
  mito_ACP_preds = predict_mito_ACPs("data/mito_ACPs.fasta", imp_ngrams, 
                                     mer_model, peptide_model),
  mito_ACP_table = get_mito_ACP_pred_table(mito_ACP_preds),
  mito_ACP_data_table = get_mito_ACP_data_table("data/mitochondrial_ACPs_table.csv"),
  property_plot = get_prop_plot(list("ACP" = pos_train_main, "AMP" = neg_train_main, "Negative" = neg_train_alt),
                                c("Hydropathy index (Kyte-Doolittle, 1982)" = "KYTJ820101", 
                                  "Net charge (Klein et al., 1984)" = "KLEP840101")),
  benchmark_table = get_benchmark_table(benchmark_peptide_preds),
  ngram_plot = get_imp_ngrams_plot(imp_ngrams_dat),
  cv_plot = get_cv_plot(cv_peptide_performance_measures),
  aa_comp_plot = get_aa_comp_plot(datasets = list("ACP" = c(pos_train_main, pos_test_main), 
                                                  "AMP" = c(neg_train_main, neg_test_main), 
                                                  "Negative" = c(neg_train_alt, neg_test_alt)))
)

make(publication_results, seed = 2938, cache = new_cache("publication_cache"))
  
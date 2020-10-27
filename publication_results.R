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

source("./functions/get_mers.R")
source("./functions/count_ngrams.R")
source("./functions/do_cv.R")
source("./functions/get_imp_ngrams.R")
source("./functions/train_models.R")
source("./functions/cdhit_data.R")
source("./functions/raw_data.R")
source("./functions/calculate_statistics.R")
source("./functions/process_data.R")
source("./functions/publication_functions.R")
source("./functions/benchmark_functions.R")


if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/Projekty/BioNgramProjects/CancerGram/publication_files/"
}
  
# This drake plan requires targets generated in analysis.R script by running 
# analysis_CancerGram drake plan. 

loadd(c(mer_df_mc, ngrams, imp_ngrams, imp_ngrams_dat, mer_model, 
        peptide_model, benchmark_mer_preds, benchmark_peptide_preds,
        pos_train_main, pos_test_main, neg_train_main, 
        neg_test_main, neg_train_alt, neg_test_alt,
        cv_mer, cv_peptide))


publication_results <- drake_plan(
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
                                                  "Negative" = c(neg_train_alt, neg_test_alt))),
  ngram_list = get_ngram_list(imp_ngrams_dat),
  gathered_acp_independent = gather_raw_data(),
  raw_acp_independent = read_raw_data("./data/all_ACPs.fasta"),
  cdhit_acp_independent = filter_cdhit(raw_acp_independent, 0.95),
  raw_amp_independent = read_amp_data(),
  cdhit_amp_independent = filter_cdhit(raw_amp_independent, 0.60),
  acp_independent = cdhit_acp_independent[which(!(cdhit_acp_independent %in% c(pos_train_main, pos_test_main, neg_train_main, neg_test_main,
                                                      neg_train_alt, neg_test_alt)))],
  amp_independent = cdhit_amp_independent[which(!(cdhit_amp_independent %in% c(pos_train_main, pos_test_main, neg_train_main, neg_test_main,
                                                      neg_train_alt, neg_test_alt)))],
  independent_dataset = write_fasta(c(acp_independent, amp_independent), "./results/independent_dataset_for_anticp_benchmark.fa"),
  independent_dataset_results = get_independent_dataset_preds(c(acp_independent, amp_independent),
                                                              "./results/independent_anticp_local_model1_0.5.csv"),
  mACPpred_benchmark_results = get_validation_mACPpred_dataset_preds(pos_mACPpred_dataset = "data/mACPpred_positive.fasta", 
                                                                     neg_mACPpred_dataset = "data/mACPpred_negative.fasta", 
                                                                     pos_validation = pos_test_main, 
                                                                     neg_validation = neg_test_main, 
                                                                     cancergram_preds = benchmark_peptide_preds, 
                                                                     mACPpred_preds = "results/mACPpred_predictions_validation.csv"),
  peptide_aa_comp = get_peptide_aa_comp(datasets = list("ACP" = c(pos_train_main, pos_test_main), 
                                                        "AMP" = c(neg_train_main, neg_test_main), 
                                                        "Negative" = c(neg_train_alt, neg_test_alt))),
  aa_comp_table_acp_amp = test_aa_comp(peptide_aa_comp, "ACP", "AMP"),
  aa_comp_table_acp_neg = test_aa_comp(peptide_aa_comp, "ACP", "Negative"),
  aa_comp_table_amp_neg = test_aa_comp(peptide_aa_comp, "AMP", "Negative"),
  anticp_benchmark_table = get_benchmark_results_table(independent_dataset_results, "AntiCP"),
  mACPpred_benchmark_table = get_benchmark_results_table(mACPpred_benchmark_results, "mACPpred")
)

make(publication_results, seed = 2938, cache = new_cache("publication_cache"))
  
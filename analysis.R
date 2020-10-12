library(dplyr)
library(magrittr)
library(drake)
library(biogram)
library(seqinr)
library(pbapply)
library(ranger)
library(cvTools)
library(visNetwork)
library(hmeasure)
library(seqR)

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/Projekty/BioNgramProjects/CancerGram/"
}

source("./functions/process_data.R")
source("./functions/get_mers.R")
source("./functions/count_ngrams.R")
source("./functions/do_cv.R")
source("./functions/get_imp_ngrams.R")
source("./functions/train_models.R")
source("./functions/calculate_statistics.R")


analysis_CancerGram <- drake_plan(pos_train_main = process_sequences("pos_train_main.txt"), # 3 sequences < 5 aa
                                  pos_test_main = process_sequences("pos_test_main.txt"), # 1 sequence < 5 aa
                                  neg_train_main = process_sequences("neg_train_main.txt"), 
                                  neg_test_main = process_sequences("neg_test_main.txt"), # 2 sequences < 5 aa
                                  neg_train_alt = process_sequences("neg_train_alternate.txt"),
                                  neg_test_alt = process_sequences("neg_test_alternate.txt"),
                                  mer_df_mc = get_target(mer_df_from_list_len_group(c(pos_train_main, neg_train_main, neg_train_alt))),
                                  ngrams = count_and_gather_ngrams(mer_df_mc,
                                                                   c(1, rep(2, 4), rep(3, 4)),
                                                                   list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
                                  cv_mer = do_cv_mc(mer_df_mc, ngrams, 0.0001),
                                  cv_peptide = do_cv_peptides_mc(cv_mer),
                                  imp_ngrams_dat = get_imp_ngrams_mc(ngrams, mer_df_mc, 0.0001),
                                  imp_ngrams = unique(unlist(unname(imp_ngrams_dat))),
                                  mer_model = train_mc_model_mers(mer_df_mc, ngrams, imp_ngrams),
                                  mer_preds = cbind(mer_df_mc, 
                                                    predict(mer_model, 
                                                            as.matrix(ngrams[, imp_ngrams]))[["predictions"]]),
                                  stats = calculate_statistics_mc(mer_preds,  c("acp", "amp", "neg")),
                                  peptide_model = train_mc_model_peptides(get_target(stats)),
                                  benchmark_mer_df = get_target(mer_df_from_list_len_group(c(pos_test_main, neg_test_main, neg_test_alt))),
                                  benchmark_ngrams = count_imp_ngrams(benchmark_mer_df, imp_ngrams),
                                  benchmark_mer_preds = cbind(benchmark_mer_df, 
                                                              predict(mer_model, 
                                                                      as.matrix(benchmark_ngrams))[["predictions"]]),
                                  benchmark_stats = get_target(calculate_statistics_mc(benchmark_mer_preds, c("acp", "amp", "neg"))),
                                  benchmark_peptide_preds = cbind(benchmark_stats[, c("source_peptide", "target")],
                                                                  predict(peptide_model, 
                                                                          benchmark_stats)[["predictions"]]))


make(analysis_CancerGram, seed = 2938)

CancerGram_model <- list(rf_mers = loadd(mer_model), rf_peptides = loadd(peptide_model), imp_features = loadd(imp_ngrams))
class(CancerGram_model) <- "cancergram_model"
save(CancerGram_model, file = "./results/CancerGram_model.rda", compress = "xz", compression_level = 9)
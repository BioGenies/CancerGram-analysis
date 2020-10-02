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

source("./functions/raw_data.R")
source("./functions/cdhit_data.R")
source("./functions/nonstandard_AMPs.R")
source("./functions/cutting_seqs.R")
source("./functions/holdouts.R")
source("./functions/writing_benchmarks.R")
source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/do_cv.R")
source("./functions/train_model_peptides.R")
source("./functions/mc_model_functions.R")
source("./functions/get_selected_features.R")
source("./functions/benchmark_functions.R")
source("./functions/degenerate_ngrams.R")
data(aaindex)

analysis_CancerGram <- drake_plan(gathered_data = gather_raw_data(),
                                  raw_data = read_raw_data("./data/all_ACPs.fasta"),
                                  cdhit_data = filter_cdhit(raw_data, 0.9),
                                  negative_data = read_and_cut(data_path, lengths(cdhit_data)),
                                  cdhit_data_ids = generate_holdout_groups(cdhit_data),
                                  negative_data_ids = generate_holdout_groups(negative_data),
                                  benchmark_file = write_benchmark(pos = cdhit_data,
                                                                   pos_id = cdhit_data_ids,
                                                                   neg = negative_data,
                                                                   neg_id = negative_data_ids),
                                  mer_df = get_mers(pos = cdhit_data,
                                                    pos_id = cdhit_data_ids,
                                                    neg = negative_data,
                                                    neg_id = negative_data_ids),
                                  binary_ngrams = count_and_gather_ngrams(mer_df,
                                                                          c(1, rep(2, 4), rep(3, 4)),
                                                                          list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
                                  cv_raw = do_cv(mer_df, binary_ngrams),
                                  pos_train_main = process_sequences("pos_train_main.txt"), # 3 sequences < 5 aa
                                  pos_test_main = process_sequences("pos_test_main.txt"), # 1 sequence < 5 aa
                                  neg_train_main = process_sequences("neg_train_main.txt"), 
                                  neg_test_main = process_sequences("neg_test_main.txt"), # 2 sequences < 5 aa
                                  neg_train_alt = process_sequences("neg_train_alternate.txt"),
                                  neg_test_alt = process_sequences("neg_test_alternate.txt"),
                                  mer_df = mutate(mer_df_from_list_len_group(c(pos_train_main, neg_train_main, neg_train_alt)),
                                                  target = case_when(grepl("pos_train_main", source_peptide) ~ "acp",
                                                                     grepl("neg_train_main", source_peptide) ~ "amp",
                                                                     grepl("neg_train_alternate", source_peptide) ~ "neg")),
                                  ngrams = count_and_gather_ngrams(mer_df,
                                                                   c(1, rep(2, 4), rep(3, 4)),
                                                                   list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
                                  imp_ngrams_dat = get_imp_ngrams_mc(ngrams, mer_df, 0.0001),
                                  imp_ngrams = unique(unlist(unname(imp_ngrams_dat))),
                                  mer_model = train_mc_model_mers(mer_df, ngrams, imp_ngrams),
                                  mer_preds = cbind(mer_df, 
                                                    predict(mer_model, 
                                                            as.matrix(ngrams[, imp_ngrams]))[["predictions"]]),
                                  stats = calculate_statistics_mc(mer_preds,  c("acp", "amp", "neg")),
                                  peptide_model = train_mc_model_peptides(get_target(stats)),
                                  benchmark_mer_df = get_target(mer_df_from_list_len_group(c(pos_test_main, neg_test_main, pos_test_alt, neg_test_alt))),
                                  benchmark_ngrams = count_imp_ngrams(benchmark_mer_df, imp_ngrams),
                                  benchmark_mer_preds = cbind(benchmark_mer_df, 
                                                              predict(mer_model, 
                                                                      as.matrix(benchmark_ngrams))[["predictions"]]),
                                  benchmark_stats = get_target(calculate_statistics_mc(benchmark_mer_preds, c("acp", "amp", "neg"))),
                                  benchmark_peptide_preds = cbind(benchmark_stats[, c("source_peptide", "target")],
                                                                  predict(peptide_model, 
                                                                          benchmark_stats)[["predictions"]]))


make(analysis_CancerGram, seed = 2938)

CancerGram_model <- list(rf_mers = mer_model, rf_peptides = peptide_model, imp_features = imp_ngrams)
class(CancerGram_model) <- "cancergram_model"
save(CancerGram_model, file = "./results/CancerGram_model.rda", compress = "xz", compression_level = 9)
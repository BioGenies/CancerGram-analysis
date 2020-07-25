library(dplyr)
library(biogram)
library(seqR)
library(ranger)
library(ggplot2)
library(tidyr)

source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/benchmark_functions.R")
source("./functions/holdouts.R")
source("./functions/mc_model_functions.R")
source("./functions/writing_benchmarks.R")
source("./functions/train_model_peptides.R")


benchmark_first_models <- drake_plan(
  # AntiCP datasets
  pos_train_main = process_sequences("pos_train_main.txt"), # 3 sequences < 5 aa
  pos_test_main = process_sequences("pos_test_main.txt"), # 1 sequence < 5 aa
  neg_train_main = process_sequences("neg_train_main.txt"), 
  neg_test_main = process_sequences("neg_test_main.txt"), # 2 sequences < 5 aa
  pos_train_alt = process_sequences("pos_train_alternate.txt"), # 4 sequences < 5 aa
  pos_test_alt = process_sequences("pos_test_alternate.txt"), 
  neg_train_alt = process_sequences("neg_train_alternate.txt"),
  neg_test_alt = process_sequences("neg_test_alternate.txt"),
  # Models on AntiCP datasets
  main_preds = train_and_test_anticp(pos_train_main, neg_train_main, pos_test_main, neg_test_main),
  alt_preds = train_and_test_anticp(pos_train_alt, neg_train_alt, pos_test_alt, neg_test_alt),
  anticp_data_pred_metrics = mapply(function(set, cutoff, name) get_metrics(set, cutoff, name), set = list(main_preds, alt_preds), 
                                    cutoff = c(rep(0.5, 2), rep(0.7, 2), rep(0.9, 2)), 
                                    name = c("main", "alt")),
  # Data for multiclass model
  acp = cdhit_data,
  amp = filter_amps(amp_full_dataset = readd(cdhit_data, path = "/home/kasia/RProjects/AmpGram-analysis/.drake"),
                    acp_dataset = acp),
  neg = negative_data,
  acp_ids = generate_holdout_groups(acp),
  amp_ids = generate_holdout_groups(amp),
  neg_ids = generate_holdout_groups(neg),
  benchmark_mc = write_benchmark_mc(acp, acp_ids, amp, amp_ids, neg, neg_ids),
  mers_mc = get_mers_mc(acp, acp_ids, amp, amp_ids, neg, neg_ids),
  ngrams_mc = count_and_gather_ngrams(mers_mc,
                                      c(1, rep(2, 4), rep(3, 4)),
                                      list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_dat_mc = get_imp_ngrams_mc(ngrams_mc, mers_mc),
  imp_ngrams_mc = unique(unlist(unname(imp_ngrams_dat_mc))),
  mer_model_mc = train_mc_model_mers(mers_mc, ngrams_mc, imp_ngrams_mc),
  mer_preds_mc = cbind(mers_mc, predict(mer_model_mc, as.matrix(ngrams_mc[, imp_ngrams_mc]))[["predictions"]]),
  stats_mc = do.call(cbind, lapply(c("acp", "amp", "neg"), function(i) 
    calculate_statistics_single(mer_preds_mc, i)))[,-c(17,18,33,34)],
  peptide_model_mc = train_mc_model_peptides(stats_mc),
  benchmark_seqs_mc = read_fasta("./results/benchmark_mc.fasta"),
  benchmark_mers_mc = mer_df_from_list_len_group(benchmark_seqs_mc),
  benchmark_ngrams_mc = count_imp_ngrams(benchmark_mers_mc, imp_ngrams_mc),
  benchmark_mer_preds_mc = cbind(benchmark_mers_mc, 
                              predict(mer_model_mc, 
                                      as.matrix(benchmark_ngrams_mc))[["predictions"]]),
  benchmark_stats_mc = do.call(cbind, lapply(c("acp", "amp", "neg"), function(i) 
    calculate_statistics_single(benchmark_mer_preds_mc, i)))[,-c(17,18,33,34)],
  benchmark_peptide_preds_mc = cbind(benchmark_stats_mc[, c("source_peptide", "target")],
                                     predict(peptide_model_mc,
                                       select(benchmark_stats_mc, -source_peptide))[["predictions"]])
  )

# saveRDS(main_preds, "./results/main_dataset_preds.RDS")
# saveRDS(alt_preds, "./results/alt_dataset_preds.RDS")


# main_preds %>% 
#   group_by(source_peptide, target) %>% 
#   summarise(frac_true = mean(pred > 0.5)) %>% 
#   ggplot(aes(x = factor(target), y = frac_true)) +
#   geom_boxplot()

# 
# ANC <- cdhit_data
# AmpGram_data <- readd(cdhit_data, path = "/home/kasia/RProjects/AmpGram-analysis/.drake")
# AMP <- AmpGram_data[which(!(AmpGram_data %in% ANC))]


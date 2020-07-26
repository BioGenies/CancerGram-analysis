library(dplyr)
library(drake)
library(biogram)
library(seqR)
library(ranger)
library(ggplot2)
library(tidyr)
library(cvTools)
library(measures)

source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/benchmark_functions.R")
source("./functions/holdouts.R")
source("./functions/mc_model_functions.R")
source("./functions/writing_benchmarks.R")
source("./functions/train_model_peptides.R")

loadd(cdhit_data)
loadd(negative_data)

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
  # Models on AntiCP datasets (mers only)
  main_preds = train_and_test_anticp(pos_train_main, neg_train_main, pos_test_main, neg_test_main),
  alt_preds = train_and_test_anticp(pos_train_alt, neg_train_alt, pos_test_alt, neg_test_alt),
  anticp_data_pred_metrics = mapply(function(set, cutoff, name) get_metrics(set, cutoff, name), set = list(main_preds, alt_preds), 
                                    cutoff = c(rep(0.5, 2), rep(0.7, 2), rep(0.9, 2)), 
                                    name = c("main", "alt")),
  # Multiclass model mers+peptides
  acp = readd(cdhit_data),
  amp = filter_amps(amp_full_dataset = readd(cdhit_data, path = "/home/kasia/RProjects/AmpGram-analysis/.drake"),
                    acp_dataset = acp),
  neg = readd(negative_data),
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
                                       select(benchmark_stats_mc, -source_peptide))[["predictions"]]),
  ### Binary models mers + peptides 
  # ACP/non-ACP our datasets
  mers_acp_neg = mutate(filter(mers_mc, target %in% c("acp", "neg")),
                        target = ifelse(target == "acp", TRUE, FALSE)),
  ngrams_acp_neg = count_and_gather_ngrams(mers_acp_neg,
                                           c(1, rep(2, 4), rep(3, 4)),
                                           list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_acp_neg = calc_imp_bigrams(mers_acp_neg, ngrams_acp_neg),
  mer_model_acp_neg = train_model_mers(mers_acp_neg, ngrams_acp_neg, imp_ngrams_acp_neg),
  mer_preds_acp_neg = mutate(mers_acp_neg,
                             pred = predict(mer_model_acp_neg, 
                                            data.frame(as.matrix(ngrams_acp_neg[, imp_ngrams_acp_neg])))[["predictions"]][, "TRUE"]),
  stats_acp_neg = calculate_statistics(mer_preds_acp_neg),
  peptide_model_acp_neg = train_model_peptides(stats_acp_neg),
  benchmark_mers_acp_neg = filter(benchmark_mers_mc, grepl("CancerPPD|AP|DRAMP|CUTTED", source_peptide)),
  benchmark_ngrams_acp_neg = count_imp_ngrams(benchmark_mers_acp_neg, imp_ngrams_acp_neg),
  benchmark_mer_preds_acp_neg = mutate(benchmark_mers_acp_neg,
                                       pred = predict(mer_model_acp_neg, 
                                                      data.frame(as.matrix(benchmark_ngrams_acp_neg[, imp_ngrams_acp_neg])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_neg = calculate_statistics(benchmark_mer_preds_acp_neg),
  benchmark_peptide_preds_acp_neg = cbind(benchmark_stats_acp_neg[, c("source_peptide", "target")],
                                          pred = predict(peptide_model_acp_neg, benchmark_stats_acp_neg)[["predictions"]][, "TRUE"]),
  
  # ACP/AMP our datasets
  mers_acp_amp = mutate(filter(mers_mc, target %in% c("acp", "amp")),
                        target = ifelse(target == "acp", TRUE, FALSE)),
  ngrams_acp_amp = count_and_gather_ngrams(mers_acp_amp,
                                           c(1, rep(2, 4), rep(3, 4)),
                                           list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_acp_amp = calc_imp_bigrams(mers_acp_amp, ngrams_acp_amp),
  mer_model_acp_amp = train_model_mers(mers_acp_amp, ngrams_acp_amp, imp_ngrams_acp_amp),
  mer_preds_acp_amp = mutate(mers_acp_amp,
                             pred = predict(mer_model_acp_amp, 
                                            data.frame(as.matrix(ngrams_acp_amp[, imp_ngrams_acp_amp])))[["predictions"]][, "TRUE"]),
  stats_acp_amp = calculate_statistics(mer_preds_acp_amp),
  peptide_model_acp_amp = train_model_peptides(stats_acp_amp),
  benchmark_mers_acp_amp = filter(benchmark_mers_mc, grepl("CancerPPD|AP|DRAMP|dbAMP", source_peptide)),
  benchmark_ngrams_acp_amp = count_imp_ngrams(benchmark_mers_acp_amp, imp_ngrams_acp_amp),
  benchmark_mer_preds_acp_amp = mutate(benchmark_mers_acp_amp,
                                       pred = predict(mer_model_acp_amp, 
                                                      data.frame(as.matrix(benchmark_ngrams_acp_amp[, imp_ngrams_acp_amp])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_amp = calculate_statistics(benchmark_mer_preds_acp_amp),
  benchmark_peptide_preds_acp_amp = cbind(benchmark_stats_acp_amp[, c("source_peptide", "target")],
                                          pred = predict(peptide_model_acp_amp, benchmark_stats_acp_amp)[["predictions"]][, "TRUE"]),
  
  # ACP/non-ACP AntiCP datasets
  mers_acp_neg_anticp = mutate(mer_df_from_list_len_group(c(pos_train_alt, neg_train_alt)),
                               target = ifelse(grepl("pos", source_peptide), TRUE, FALSE)),
  ngrams_acp_neg_anticp = count_and_gather_ngrams(mers_acp_neg_anticp,
                                                  c(1, rep(2, 4), rep(3, 4)),
                                                  list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_acp_neg_anticp = calc_imp_bigrams(mers_acp_neg_anticp, ngrams_acp_neg_anticp),
  mer_model_acp_neg_anticp = train_model_mers(mers_acp_neg_anticp, ngrams_acp_neg_anticp, imp_ngrams_acp_neg_anticp),
  mer_preds_acp_neg_anticp = mutate(mers_acp_neg_anticp,
                                    pred = predict(mer_model_acp_neg_anticp, 
                                                   data.frame(as.matrix(ngrams_acp_neg_anticp[, imp_ngrams_acp_neg_anticp])))[["predictions"]][, "TRUE"]),
  stats_acp_neg_anticp = calculate_statistics(mer_preds_acp_neg_anticp),
  peptide_model_acp_neg_anticp = train_model_peptides(stats_acp_neg_anticp),
  benchmark_mers_acp_neg_anticp = mutate(mer_df_from_list_len_group(c(pos_test_alt, neg_test_alt)),
                                         target = ifelse(grepl("pos", source_peptide), TRUE, FALSE)),
  benchmark_ngrams_acp_neg_anticp = count_imp_ngrams(benchmark_mers_acp_neg_anticp, imp_ngrams_acp_neg_anticp),
  benchmark_mer_preds_acp_neg_anticp = mutate(benchmark_mers_acp_neg_anticp,
                                              pred = predict(mer_model_acp_neg_anticp, 
                                                             data.frame(as.matrix(benchmark_ngrams_acp_neg_anticp[, imp_ngrams_acp_neg_anticp])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_neg_anticp = calculate_statistics(benchmark_mer_preds_acp_neg_anticp),
  benchmark_peptide_preds_acp_neg_anticp = cbind(benchmark_stats_acp_neg_anticp[, c("source_peptide", "target")],
                                                 pred = predict(peptide_model_acp_neg_anticp, benchmark_stats_acp_neg_anticp)[["predictions"]][, "TRUE"]),
  
  # ACP/AMP AntiCP datasets
  mers_acp_amp_anticp = mutate(mer_df_from_list_len_group(c(pos_train_main, neg_train_main)),
                               target = ifelse(grepl("pos", source_peptide), TRUE, FALSE)),
  ngrams_acp_amp_anticp = count_and_gather_ngrams(mers_acp_amp_anticp,
                                           c(1, rep(2, 4), rep(3, 4)),
                                           list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_acp_amp_anticp = calc_imp_bigrams(mers_acp_amp_anticp, ngrams_acp_amp_anticp),
  mer_model_acp_amp_anticp = train_model_mers(mers_acp_amp_anticp, ngrams_acp_amp_anticp, imp_ngrams_acp_amp_anticp),
  mer_preds_acp_amp_anticp = mutate(mers_acp_amp_anticp,
                             pred = predict(mer_model_acp_amp_anticp, 
                                            data.frame(as.matrix(ngrams_acp_amp_anticp[, imp_ngrams_acp_amp_anticp])))[["predictions"]][, "TRUE"]),
  stats_acp_amp_anticp = calculate_statistics(mer_preds_acp_amp_anticp),
  peptide_model_acp_amp_anticp = train_model_peptides(stats_acp_amp_anticp),
  benchmark_mers_acp_amp_anticp = mutate(mer_df_from_list_len_group(c(pos_test_main, neg_test_main)),
                                  target = ifelse(grepl("pos", source_peptide), TRUE, FALSE)),
  benchmark_ngrams_acp_amp_anticp = count_imp_ngrams(benchmark_mers_acp_amp_anticp, imp_ngrams_acp_amp_anticp),
  benchmark_mer_preds_acp_amp_anticp = mutate(benchmark_mers_acp_amp_anticp,
                                       pred = predict(mer_model_acp_amp_anticp, 
                                                      data.frame(as.matrix(benchmark_ngrams_acp_amp_anticp[, imp_ngrams_acp_amp_anticp])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_amp_anticp = calculate_statistics(benchmark_mer_preds_acp_amp_anticp),
  benchmark_peptide_preds_acp_amp_anticp = cbind(benchmark_stats_acp_amp_anticp[, c("source_peptide", "target")],
                                          pred = predict(peptide_model_acp_amp_anticp, benchmark_stats_acp_amp_anticp)[["predictions"]][, "TRUE"])
  )

make(benchmark_first_models, seed = 2938)

# saveRDS(main_preds, "./results/main_dataset_preds.RDS")
# saveRDS(alt_preds, "./results/alt_dataset_preds.RDS")


# main_preds %>% 
#   group_by(source_peptide, target) %>% 
#   summarise(frac_true = mean(pred > 0.5)) %>% 
#   ggplot(aes(x = factor(target), y = frac_true)) +
#   geom_boxplot()

### Performance measures

mc <- readd(benchmark_peptide_preds_mc)
acp_neg <- readd(benchmark_peptide_preds_acp_neg)
acp_amp <- readd(benchmark_peptide_preds_acp_amp)
acp_neg_anticp <- readd(benchmark_peptide_preds_acp_neg_anticp)
acp_amp_anticp <- readd(benchmark_peptide_preds_acp_amp_anticp)

res_list <- list("acp_neg" = acp_neg, "acp_amp" = acp_amp, "acp_neg_anticp" = acp_neg_anticp, "acp_amp_anticp" = acp_amp_anticp)

all_res <- lapply(1:length(res_list), function(i) {
  mutate(res_list[[i]], model = names(res_list)[i])
}) %>% bind_rows()


mutate(all_res, target = factor(target),
       decision = factor(ifelse(pred > 0.5, TRUE, FALSE))) %>% 
  group_by(model) %>% 
  summarise(MCC = mlr3measures::mcc(target, decision, "TRUE"),
           Precision = mlr3measures::precision(target, decision, "TRUE"),
           Sensitivity = mlr3measures::sensitivity(target, decision, "TRUE"),
           Specificity = mlr3measures::specificity(target, decision, "TRUE"),
           Accuracy = mlr3measures::acc(target, decision),
           AUC = mlr3measures::auc(target, pred, "TRUE")) 

# Multiclass model
multiclass.AUNP(mc[, c("acp", "amp", "neg")], mc[["target"]])
multiclass.AU1U(mc[, c("acp", "amp", "neg")], mc[["target"]])
multiclass.AUNU(mc[, c("acp", "amp", "neg")], mc[["target"]])

mc_decisions <- mutate(mc, decision = case_when(acp > amp & acp > neg  ~ "acp",
                                                amp > acp & amp > neg ~ "amp",
                                                neg > amp & neg > acp ~ "neg"))

ACC(mc_decisions[["target"]], mc_decisions[["decision"]])
table(mc_decisions[,c(2,6)])
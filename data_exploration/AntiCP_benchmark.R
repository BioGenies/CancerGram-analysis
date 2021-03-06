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

if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/Projekty/BioNgramProjects/CancerGram/"
}

source("./functions/raw_data.R")
source("./functions/cdhit_data.R")
source("./functions/cutting_seqs.R")
source("./functions/holdouts.R")
source("./functions/writing_benchmarks.R")
source("../functions/get_mers.R")
source("./functions/do_cv.R")
source("../functions/count_ngrams.R")
source("./functions/holdouts.R")
source("./functions/exploration_functions")
source("./functions/writing_benchmarks.R")
source("../functions/calculate_statistics.R")
source("../functions/do_cv.R")
source("../functions/get_imp_ngrams.R")
source("../functions/process_data.R")
source("../functions/train_models.R")

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
  
  # Multiclass model mers+peptides our datasets
  gathered_acp_data = gather_raw_data(),
  raw_acp_data = read_raw_data("../data/all_ACPs.fasta"),
  cdhit_acp_data = filter_cdhit(raw_acp_data, 0.95),
  raw_amp_data = read_amp_data(),
  cdhit_amp_data = filter_cdhit(raw_amp_data, 0.95),
  amp_filtered_data = filter_amps(cdhit_amp_data, cdhit_acp_data),
  neg_data = read_and_cut(data_path, lengths(cdhit_acp_data)),
  cdhit_acp_data_ids = generate_holdout_groups(cdhit_acp_data),
  amp_filtered_data_ids = generate_holdout_groups(amp_filtered_data),
  neg_data_ids = generate_holdout_groups(neg_data),
  benchmark_file = write_benchmark_mc(acp = cdhit_acp_data,
                                      acp_id = cdhit_acp_data_ids,
                                      amp = amp_filtered_data,
                                      amp_id = amp_filtered_data_ids,
                                      neg = neg_data,
                                      neg_id = neg_data_ids),
  mer_df_mc = get_mers_mc(acp = cdhit_acp_data,
                          acp_id = cdhit_acp_data_ids,
                          amp = amp_filtered_data,
                          amp_id = amp_filtered_data_ids,
                          neg = neg_data,
                          neg_id = neg_data_ids),
  ngrams_mc = count_and_gather_ngrams(mer_df_mc,
                                      c(1, rep(2, 4), rep(3, 4)),
                                      list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  cv_raw_mc = do_cv_mc(mer_df_mc, ngrams_mc, 0.05),
  benchmark_mc = write_benchmark_mc(acp = cdhit_acp_data,
                                    acp_id = cdhit_acp_data_ids,
                                    amp = amp_filtered_data,
                                    amp_id = amp_filtered_data_ids,
                                    neg = neg_data,
                                    neg_id = neg_data_ids),
  imp_ngrams_dat_mc = get_imp_ngrams_mc(ngrams_mc, mer_df_mc, 0.05),
  imp_ngrams_mc = unique(unlist(unname(imp_ngrams_dat_mc))),
  mer_model_mc = train_mc_model_mers(mer_df_mc, ngrams_mc, imp_ngrams_mc),
  mer_preds_mc = cbind(mer_df_mc, predict(mer_model_mc, as.matrix(ngrams_mc[, imp_ngrams_mc]))[["predictions"]]),
  stats_mc = mutate(calculate_statistics_mc(mer_preds_mc, c("acp", "amp", "neg")),
                    target = factor(case_when(grepl("CUTTED", source_peptide) ~ "neg",
                                              grepl("dbAMP", source_peptide) ~ "amp",
                                              grepl("CancerPPD|AP|DRAMP", source_peptide) ~ "acp"))),
  peptide_model_mc = train_mc_model_peptides(stats_mc),
  benchmark_seqs_mc = read_fasta("./results/benchmark_mc.fasta"),
  benchmark_mers_mc = mer_df_from_list_len_group(benchmark_seqs_mc),
  benchmark_ngrams_mc = count_imp_ngrams(benchmark_mers_mc, imp_ngrams_mc),
  benchmark_mer_preds_mc = cbind(benchmark_mers_mc, 
                                 predict(mer_model_mc, 
                                         as.matrix(benchmark_ngrams_mc))[["predictions"]]),
  benchmark_stats_mc = mutate(calculate_statistics_mc(benchmark_mer_preds_mc, c("acp", "amp", "neg")),
                              target = factor(case_when(grepl("CUTTED", source_peptide) ~ "neg",
                                                        grepl("dbAMP", source_peptide) ~ "amp",
                                                        grepl("CancerPPD|AP|DRAMP", source_peptide) ~ "acp"))),
  benchmark_peptide_preds_mc = cbind(benchmark_stats_mc[, c("source_peptide", "target")],
                                     predict(peptide_model_mc, benchmark_stats_mc)[["predictions"]]),
  
  # Multiclass model on AntiCP datasets
  mers_mc_anticp = mutate(mer_df_from_list_len_group(c(pos_train_main, neg_train_main, neg_train_alt)),
                          target = case_when(grepl("pos_train_main", source_peptide) ~ "acp",
                                             grepl("neg_train_main", source_peptide) ~ "amp",
                                             grepl("neg_train_alternate", source_peptide) ~ "neg")),
  ngrams_mc_anticp = count_and_gather_ngrams(mers_mc_anticp,
                                             c(1, rep(2, 4), rep(3, 4)),
                                             list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_dat_mc_anticp = get_imp_ngrams_mc(ngrams_mc_anticp, mers_mc_anticp),
  imp_ngrams_mc_anticp = unique(unlist(unname(imp_ngrams_dat_mc_anticp))),
  mer_model_mc_anticp = train_mc_model_mers(mers_mc_anticp, ngrams_mc_anticp, imp_ngrams_mc_anticp),
  mer_preds_mc_anticp = cbind(mers_mc_anticp, predict(mer_model_mc_anticp, as.matrix(ngrams_mc_anticp[, imp_ngrams_mc_anticp]))[["predictions"]]),
  stats_mc_anticp = calculate_statistics_mc(mer_preds_mc_anticp,  c("acp", "amp", "neg")),
  peptide_model_mc_anticp = train_mc_model_peptides(mutate(stats_mc_anticp,
                                                           target = factor(case_when(grepl("pos_train_main", source_peptide) ~ "acp",
                                                                                     grepl("neg_train_main", source_peptide) ~ "amp",
                                                                                     grepl("neg_train_alternate", source_peptide) ~ "neg")))),
  benchmark_mers_mc_anticp = mutate(mer_df_from_list_len_group(c(pos_test_main, neg_test_main, pos_test_alt, neg_test_alt)),
                                    target = factor(case_when(grepl("pos_test", source_peptide) ~ "acp",
                                                              grepl("neg_test_main", source_peptide) ~ "amp",
                                                              grepl("neg_test_alternate", source_peptide) ~ "neg"))),
  benchmark_ngrams_mc_anticp = count_imp_ngrams(benchmark_mers_mc_anticp, imp_ngrams_mc_anticp),
  benchmark_mer_preds_mc_anticp = cbind(benchmark_mers_mc_anticp, 
                                        predict(mer_model_mc_anticp, 
                                                as.matrix(benchmark_ngrams_mc_anticp))[["predictions"]]),
  benchmark_stats_mc_anticp = mutate(calculate_statistics_mc(benchmark_mer_preds_mc_anticp, c("acp", "amp", "neg")),
                                     target = factor(case_when(grepl("pos_test", source_peptide) ~ "acp",
                                                               grepl("neg_test_main", source_peptide) ~ "amp",
                                                               grepl("neg_test_alternate", source_peptide) ~ "neg"))),
  benchmark_peptide_preds_mc_anticp = cbind(benchmark_stats_mc_anticp[, c("source_peptide", "target")],
                                            predict(peptide_model_mc_anticp, benchmark_stats_mc_anticp)[["predictions"]]),
  ### Binary models mers + peptides 
  # ACP/non-ACP our datasets
  mers_acp_neg = mutate(filter(mer_df_mc, target %in% c("acp", "neg")),
                        target = ifelse(target == "acp", TRUE, FALSE)),
  ngrams_acp_neg = count_and_gather_ngrams(mers_acp_neg,
                                           c(1, rep(2, 4), rep(3, 4)),
                                           list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_acp_neg = calc_imp_bigrams(mers_acp_neg, ngrams_acp_neg),
  mer_model_acp_neg = train_model_mers(mers_acp_neg, ngrams_acp_neg, imp_ngrams_acp_neg),
  mer_preds_acp_neg = mutate(mers_acp_neg,
                             pred = predict(mer_model_acp_neg, 
                                            data.frame(as.matrix(ngrams_acp_neg[, imp_ngrams_acp_neg])))[["predictions"]][, "TRUE"]),
  stats_acp_neg = mutate(calculate_statistics(mer_preds_acp_neg),
                         target = factor(ifelse(grepl("CUTTED", source_peptide), FALSE, TRUE))),
  peptide_model_acp_neg = train_model_peptides(stats_acp_neg),
  benchmark_mers_acp_neg = filter(benchmark_mers_mc, grepl("CancerPPD|AP|DRAMP|CUTTED", source_peptide)),
  benchmark_ngrams_acp_neg = count_imp_ngrams(benchmark_mers_acp_neg, imp_ngrams_acp_neg),
  benchmark_mer_preds_acp_neg = mutate(benchmark_mers_acp_neg,
                                       pred = predict(mer_model_acp_neg, 
                                                      data.frame(as.matrix(benchmark_ngrams_acp_neg[, imp_ngrams_acp_neg])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_neg = mutate(calculate_statistics(benchmark_mer_preds_acp_neg),
                                   target = factor(ifelse(grepl("CUTTED", source_peptide), FALSE, TRUE))),
  benchmark_peptide_preds_acp_neg = mutate(benchmark_stats_acp_neg[, c("source_peptide", "target")],
                                           pred = predict(peptide_model_acp_neg, benchmark_stats_acp_neg)[["predictions"]][, "TRUE"]),
  
  # ACP/AMP our datasets
  mers_acp_amp = mutate(filter(mer_df_mc, target %in% c("acp", "amp")),
                        target = ifelse(target == "acp", TRUE, FALSE)),
  ngrams_acp_amp = count_and_gather_ngrams(mers_acp_amp,
                                           c(1, rep(2, 4), rep(3, 4)),
                                           list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
  imp_ngrams_acp_amp = calc_imp_bigrams(mers_acp_amp, ngrams_acp_amp),
  mer_model_acp_amp = train_model_mers(mers_acp_amp, ngrams_acp_amp, imp_ngrams_acp_amp),
  mer_preds_acp_amp = mutate(mers_acp_amp,
                             pred = predict(mer_model_acp_amp, 
                                            data.frame(as.matrix(ngrams_acp_amp[, imp_ngrams_acp_amp])))[["predictions"]][, "TRUE"]),
  stats_acp_amp = mutate(calculate_statistics(mer_preds_acp_amp),
                         target = factor(ifelse(grepl("dbAMP", source_peptide), FALSE, TRUE))),
  peptide_model_acp_amp = train_model_peptides(stats_acp_amp),
  benchmark_mers_acp_amp = filter(benchmark_mers_mc, grepl("CancerPPD|AP|DRAMP|dbAMP", source_peptide)),
  benchmark_ngrams_acp_amp = count_imp_ngrams(benchmark_mers_acp_amp, imp_ngrams_acp_amp),
  benchmark_mer_preds_acp_amp = mutate(benchmark_mers_acp_amp,
                                       pred = predict(mer_model_acp_amp, 
                                                      data.frame(as.matrix(benchmark_ngrams_acp_amp[, imp_ngrams_acp_amp])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_amp = mutate(calculate_statistics(benchmark_mer_preds_acp_amp),
                                   target = factor(ifelse(grepl("dbAMP", source_peptide), FALSE, TRUE))),
  benchmark_peptide_preds_acp_amp = mutate(benchmark_stats_acp_amp[, c("source_peptide", "target")],
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
  stats_acp_neg_anticp = mutate(calculate_statistics(mer_preds_acp_neg_anticp),
                                target = factor(ifelse(grepl("pos", source_peptide), TRUE, FALSE))),
  peptide_model_acp_neg_anticp = train_model_peptides(stats_acp_neg_anticp),
  benchmark_mers_acp_neg_anticp = mutate(mer_df_from_list_len_group(c(pos_test_alt, neg_test_alt)),
                                         target = factor(ifelse(grepl("pos", source_peptide), TRUE, FALSE))),
  benchmark_ngrams_acp_neg_anticp = count_imp_ngrams(benchmark_mers_acp_neg_anticp, imp_ngrams_acp_neg_anticp),
  benchmark_mer_preds_acp_neg_anticp = mutate(benchmark_mers_acp_neg_anticp,
                                              pred = predict(mer_model_acp_neg_anticp, 
                                                             data.frame(as.matrix(benchmark_ngrams_acp_neg_anticp[, imp_ngrams_acp_neg_anticp])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_neg_anticp = mutate(calculate_statistics(benchmark_mer_preds_acp_neg_anticp),
                                          target = factor(ifelse(grepl("pos", source_peptide), TRUE, FALSE))),
  benchmark_peptide_preds_acp_neg_anticp = mutate(benchmark_stats_acp_neg_anticp[, c("source_peptide", "target")],
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
  stats_acp_amp_anticp = mutate(calculate_statistics(mer_preds_acp_amp_anticp),
                                target = factor(ifelse(grepl("pos", source_peptide), TRUE, FALSE))),
  peptide_model_acp_amp_anticp = train_model_peptides(stats_acp_amp_anticp),
  benchmark_mers_acp_amp_anticp = mutate(mer_df_from_list_len_group(c(pos_test_main, neg_test_main)),
                                         target = factor(ifelse(grepl("pos", source_peptide), TRUE, FALSE))),
  benchmark_ngrams_acp_amp_anticp = count_imp_ngrams(benchmark_mers_acp_amp_anticp, imp_ngrams_acp_amp_anticp),
  benchmark_mer_preds_acp_amp_anticp = mutate(benchmark_mers_acp_amp_anticp,
                                              pred = predict(mer_model_acp_amp_anticp, 
                                                             data.frame(as.matrix(benchmark_ngrams_acp_amp_anticp[, imp_ngrams_acp_amp_anticp])))[["predictions"]][, "TRUE"]),
  benchmark_stats_acp_amp_anticp = mutate(calculate_statistics(benchmark_mer_preds_acp_amp_anticp),
                                          target = factor(ifelse(grepl("pos", source_peptide), TRUE, FALSE))),
  benchmark_peptide_preds_acp_amp_anticp = mutate(benchmark_stats_acp_amp_anticp[, c("source_peptide", "target")],
                                                  pred = predict(peptide_model_acp_amp_anticp, benchmark_stats_acp_amp_anticp)[["predictions"]][, "TRUE"]),
  
  # Create validation dataset composed of both our and AntiCP benchmark datasets to compare multiclass models:
  validation_dataset = get_validation_dataset(),
  validation_mer_df = mutate(mer_df_from_list_len_group(validation_dataset),
                             target = case_when(grepl("pos_test_main|Cancer|AP|DRAMP", source_peptide) ~ "acp",
                                                grepl("neg_test_main|dbAMP", source_peptide) ~ "amp",
                                                grepl("neg_test_alternate|CUTTED", source_peptide) ~ "neg")),
  validation_anticp_ngrams = count_imp_ngrams(validation_mer_df, imp_ngrams_mc_anticp),
  validation_anticp_mer_preds = cbind(validation_mer_df, 
                                   predict(mer_model_mc_anticp, 
                                           as.matrix(validation_anticp_ngrams))[["predictions"]]),
  validation_anticp_stats = mutate(calculate_statistics_mc(validation_anticp_mer_preds,  c("acp", "amp", "neg")),
                                   target = case_when(grepl("pos_test_main|Cancer|AP|DRAMP", source_peptide) ~ "acp",
                                                      grepl("neg_test_main|dbAMP", source_peptide) ~ "amp",
                                                      grepl("neg_test_alternate|CUTTED", source_peptide) ~ "neg")),
  validation_anticp_peptide_preds = cbind(validation_anticp_stats[, c("source_peptide", "target")],
                                            predict(peptide_model_mc_anticp, validation_anticp_stats)[["predictions"]]),
  validation_our_ngrams = count_imp_ngrams(validation_mer_df, imp_ngrams_mc),
  validation_our_mer_preds = cbind(validation_mer_df, 
                                      predict(mer_model_mc, 
                                              as.matrix(validation_our_ngrams))[["predictions"]]),
  validation_our_stats = mutate(calculate_statistics_mc(validation_our_mer_preds,  c("acp", "amp", "neg")),
                                   target = case_when(grepl("pos_test_main|Cancer|AP|DRAMP", source_peptide) ~ "acp",
                                                      grepl("neg_test_main|dbAMP", source_peptide) ~ "amp",
                                                      grepl("neg_test_alternate|CUTTED", source_peptide) ~ "neg")),
  validation_our_peptide_preds = cbind(validation_our_stats[, c("source_peptide", "target")],
                                          predict(peptide_model_mc, validation_our_stats)[["predictions"]])
  
)

make(benchmark_first_models, seed = 2938, cache = new_cache("benchmark_cache"))

# saveRDS(main_preds, "./results/main_dataset_preds.RDS")
# saveRDS(alt_preds, "./results/alt_dataset_preds.RDS")


# main_preds %>% 
#   group_by(source_peptide, target) %>% 
#   summarise(frac_true = mean(pred > 0.5)) %>% 
#   ggplot(aes(x = factor(target), y = frac_true)) +
#   geom_boxplot()

### Performance measures

mc <- readd(benchmark_peptide_preds_mc)
mc_anticp <- readd(benchmark_peptide_preds_mc_anticp)
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
            Sensitivity = mlr3measures::sensitivity(target, decision, "TRUE")*100,
            Specificity = mlr3measures::specificity(target, decision, "TRUE")*100,
            Accuracy = mlr3measures::acc(target, decision)*100,
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



all_res_mc <- bind_rows(mutate(mc, model = "mc"),
                        mutate(mc_anticp, model = "mc_anticp")) %>% 
  mutate(decision = factor(ifelse(acp > amp & acp > neg, "TRUE", "FALSE")))

# Porównanie ACP/AMP
filter(all_res_mc, grepl("AP|Cancer|DRAMP|dbAMP|test_main", source_peptide)) %>% 
  mutate(target = factor(ifelse(grepl("AP|Cancer|DRAMP|pos_", source_peptide), TRUE, FALSE))) %>% 
  group_by(model) %>% 
  summarise(MCC = mlr3measures::mcc(target, decision, "TRUE"),
            Sensitivity = mlr3measures::sensitivity(target, decision, "TRUE")*100,
            Specificity = mlr3measures::specificity(target, decision, "TRUE")*100,
            Accuracy = mlr3measures::acc(target, decision)*100)

# Porównanie ACP/neg
filter(all_res_mc, grepl("AP|Cancer|DRAMP|CUTTED|pos_test|test_alt", source_peptide)) %>% 
  mutate(target = factor(ifelse(grepl("AP|Cancer|DRAMP|pos_", source_peptide), TRUE, FALSE))) %>% 
  group_by(model) %>% 
  summarise(MCC = mlr3measures::mcc(target, decision, "TRUE"),
            Sensitivity = mlr3measures::sensitivity(target, decision, "TRUE")*100,
            Specificity = mlr3measures::specificity(target, decision, "TRUE")*100,
            Accuracy = mlr3measures::acc(target, decision)*100)


### Find improperly predicted peptides (MC model on AntiCP data)
readd(benchmark_peptide_preds_mc_anticp) %>% 
  mutate(decision = case_when(acp > amp & acp > neg  ~ "acp",
                              amp > acp & amp > neg ~ "amp",
                              neg > amp & neg > acp ~ "neg")) %>% 
  filter(decision != target) %>% 
  write.csv("improperly_predicted.csv", row.names = FALSE)

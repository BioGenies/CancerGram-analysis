get_decision_mc <- function(preds) {
  mutate(preds, decision = case_when(acp > amp & acp > neg  ~ "acp",
                                     amp > acp & amp > neg ~ "amp",
                                     neg > amp & neg > acp ~ "neg"))
}

calc_cv_performance <- function(preds) {
  res <- get_decision_mc(preds)
  lapply(unique(res[["fold"]]), function(ith_fold) {
    dat <- filter(res, fold == ith_fold)
    data.frame(
      fold = ith_fold,
      Accuracy = ACC(dat[["target"]], dat[["decision"]]),
      AU1U = multiclass.AU1U(dat[, c("acp", "amp", "neg")], dat[["target"]]),
      KapS = KAPPA(dat[["target"]], dat[["decision"]]),
      stringsAsFactors = FALSE)
  }) %>% bind_rows()
}



get_cv_pred_table <- function(cv_mer_performance, cv_peptide_performance) {
  tab <- mutate(cv_mer_performance,
                Layer = "Mer layer") %>% 
    bind_rows(mutate(cv_peptide_performance,
                     Layer = "Peptide layer")) %>% 
    pivot_longer(Accuracy:KapS, names_to = "Measure", values_to = "Value") %>% 
    group_by(Layer, Measure) %>%
    summarise(mean = mean(Value),
              sd = sd(Value)) %>% 
    mutate(value = paste0(round(mean, 2), " (+/-", round(sd, 3), ")")) %>% 
    select(-c(mean, sd)) %>%
    pivot_wider(names_from = Layer, values_from = value)
  write.csv(tab, paste0(data_path, "cv_results.csv"), row.names = FALSE)
  xtable(tab, caption = "", label = "", align = "cccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
      writeLines(paste0(data_path, "cv_results.txt"))
}

get_datasets_table <- function(acp_train, acp_test, amp_train, amp_test, neg_train, neg_test) {
  df <- data.frame(Dataset = c("Train-test", "Benchmark"),
             ACP = c(length(acp_train), length(acp_test)),
             AMP = c(length(amp_train), length(amp_test)),
             Negative = c(length(neg_train), length(neg_test)))
  write.csv(df, paste0(data_path, "datasets_table.csv"), row.names = FALSE)
  xtable(df, caption = "", label = "", align = "ccccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "datasets_table.txt"))
}


predict_mito_ACPs <- function(mito_acp_file, imp_ngrams, mer_model, peptide_model) {
  peptides <- read_fasta(mito_acp_file)
  mers <- mer_df_from_list(peptides)
  ngrams <- count_imp_ngrams(mers, imp_ngrams)
  mer_preds <- cbind(mers, predict(mer_model, as.matrix(ngrams[, imp_ngrams]))[["predictions"]])
  stats <- calculate_statistics_mc(mer_preds,  c("acp", "amp", "neg"))
  peptide_preds <- cbind(stats[, c("source_peptide")],
                         predict(peptide_model, stats)[["predictions"]])
  return(list("mer_predictions" = mer_preds, "peptide_predictions" = peptide_preds))
}


get_mito_ACP_pred_table <- function(mito_ACP_preds) {
  res <- mito_ACP_preds[["peptide_predictions"]] %>% 
    data.frame(stringsAsFactors=FALSE) %>%
    mutate(Decision = case_when(acp > amp & acp > neg  ~ "ACP",
                                       amp > acp & amp > neg ~ "AMP",
                                       neg > amp & neg > acp ~ "Negative")) %>% 
    mutate_at(vars(acp, amp, neg), ~round(as.numeric(.), 3))
  colnames(res) <- c("Peptide", "ACP", "AMP", "Negative", "Decision")
  write.csv(res, paste0(data_path, "mito_ACP_preds.csv"), row.names = FALSE)
  xtable(res, caption = "", label = "Tab:mito_ACP_preds", align = "cccccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "mito_ACP_preds.txt"))
}

encode_seq <- function(x, property) {
  sapply(x, function(ith_seq) {
    mean(aaprop[property, tolower(ith_seq)])
  })
}

calc_prop_values <- function(seq_list, prop_vec) {
  lapply(names(prop_vec), function(ith_prop_type) {
    lapply(names(seq_list), function(ith_seq_type) {
      data.frame(seq_type = ith_seq_type, prop_type = ith_prop_type, 
                 value = encode_seq(seq_list[[ith_seq_type]], prop_vec[ith_prop_type]),
                 stringsAsFactors = FALSE)
    }) %>% bind_rows
  }) %>% bind_rows()
}

get_prop_plot <- function(seq_list, prop_vec) {
  p <- calc_prop_values(seq_list, prop_vec) %>% 
    group_by(prop_type) %>% 
    mutate(row = row_number()) %>% 
    pivot_wider(names_from = prop_type, values_from = value) %>% 
    ggplot(aes(x = `Hydropathy index (Kyte-Doolittle, 1982)`, y = `Net charge (Klein et al., 1984)`, fill = seq_type)) +
    stat_density_2d(aes(alpha = ..level..), geom = "polygon", color = "black", size = 0.4) +
    scale_fill_manual("Dataset", values = c("#ed463d", "#ffc745", "#c3dae8")) +
    facet_wrap(~seq_type) +
    guides(alpha = FALSE) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(plot = p, filename = paste0(data_path, "prop_plot.eps"), device = cairo_ps, height = 5, width = 8)
  p
}

get_mito_ACP_data_table <- function(filename) {
  tab <- read.csv(filename, stringsAsFactors = FALSE) %>% 
    arrange(Peptide)
  xtab <- mutate(tab, Reference = paste0("citep{", Reference, "}"))
  write.csv(tab, paste0(data_path, "mito_ACP_data_table.csv"), row.names = FALSE)
  xtable(xtab, caption = "", label = "", align = "cccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "mito_ACP_data_table.txt"))
}


get_benchmark_table <- function(benchmark_res) {
  dat <- get_decision_mc(benchmark_res)
  res <- data.frame(Accuracy = ACC(dat[["target"]], dat[["decision"]]),
                    AU1U = multiclass.AU1U(dat[, c("acp", "amp", "neg")], dat[["target"]]),
                    KapS = KAPPA(dat[["target"]], dat[["decision"]]),
                    stringsAsFactors = FALSE) %>% 
    pivot_longer(Accuracy:KapS, names_to = "Measure", values_to = "Value")
  write.csv(res, paste0(data_path, "benchmark_res.csv"), row.names = FALSE)
  xtable(res, caption = "", label = "Tab:benchmark", align = "ccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "benchmark_res.txt"))
}


get_imp_ngrams_plot <- function(imp_ngrams_dat) {
  res <- lapply(names(imp_ngrams_dat), function(ith_comb) {
    all <- paste0(imp_ngrams_dat[[ith_comb]], collapse = "")
    aa <- gsub("1|2|3|0|_", "", all)
    aa <- gsub(".", "", aa, fixed = TRUE)
    data.frame(table(strsplit(aa, "")), stringsAsFactors = FALSE) %>% 
      mutate(Freq = Freq/nchar(aa),
             comb = ith_comb)
  }) %>% bind_rows() 
  colnames(res) <- c("Amino acid", "Frequency", "Comparison")
  
  p <- ggplot(res, aes(x = `Amino acid`, y = Frequency, fill = Comparison)) +
    geom_col(position = "dodge") +
    scale_fill_manual("Comparison", values = c("#f0a87a", "#c7bae1", "#b9dfc8"), 
                      labels = c("ACP/AMP", "ACP/Negative", "AMP/Negative")) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(plot = p, filename = paste0(data_path, "ngram_plot.eps"), device = cairo_ps, height = 5, width = 9)
}


get_cv_plot <- function(cv_perf) {
  p <- pivot_longer(cv_perf, Accuracy:KapS, names_to = "Measure", values_to = "Value") %>% 
    ggplot(aes(x = Measure, y = Value)) +
    geom_point(position = "jitter") +
    facet_wrap(~Measure, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  ggsave(plot = p, filename = paste0(data_path, "cv_plot.eps"), device = cairo_ps, height = 4, width = 6)
}


get_aa_comp_plot <- function(datasets) {
  
  res <- lapply(names(datasets), function(ith_dataset) {
    aac <- data.frame(table(unlist(datasets[ith_dataset]))) %>% 
      setNames(c("Amino acid", "Frequency")) 
    mutate(aac, 
           Frequency = Frequency/sum(Frequency),
           Dataset = ith_dataset)
  }) %>% bind_rows() 
  
  p <- ggplot(res, aes(x = `Amino acid`, y = Frequency, fill = Dataset)) +
    geom_col(position = "dodge") +
    scale_fill_manual("Dataset", values = c("#ed463d", "#ffc745", "#c3dae8")) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(plot = p, filename = paste0(data_path, "aa_comp_plot.eps"), device = cairo_ps, height = 5, width = 9)
}


get_peptide_aa_comp <- function(datasets) {
  comp <- lapply(names(datasets), function(ith_dataset) {
    lapply(names(datasets[[ith_dataset]]), function(ith_prot) {
      aac <- data.frame(table(factor(unlist(datasets[[ith_dataset]][[ith_prot]]), 
                                     levels = c("A", "C", "D", "E", "F", "G", "H", 
                                                "I", "K", "L", "M", "N", "P", "Q", 
                                                "R", "S", "T", "V", "W", "Y")))) %>% 
        setNames(c("Amino acid", "Frequency")) 
      mutate(aac, 
             Frequency = Frequency/sum(Frequency),
             Dataset = ith_dataset)
    }) %>% bind_rows()
  }) %>% bind_rows() 
}

test_aa_comp <- function(peptide_aa_comp, dataset1, dataset2) {
  test_res <- lapply(unique(peptide_aa_comp[["Amino acid"]]), function(ith_aa) {
    data.frame(aa = ith_aa,
               pval = wilcox.test(x = filter(peptide_aa_comp, `Amino acid` == ith_aa, Dataset == dataset1)[["Frequency"]],
                                  y = filter(peptide_aa_comp, `Amino acid` == ith_aa, Dataset == dataset2)[["Frequency"]],
                                  correct = FALSE)[["p.value"]],
               mean1 = round(mean(filter(peptide_aa_comp, `Amino acid` == ith_aa, Dataset == dataset1)[["Frequency"]])*100, 2),
               mean2 = round(mean(filter(peptide_aa_comp, `Amino acid` == ith_aa, Dataset == dataset2)[["Frequency"]])*100, 2),
               median1 = round(median(filter(peptide_aa_comp, `Amino acid` == ith_aa, Dataset == dataset1)[["Frequency"]])*100, 2),
               median2 = round(median(filter(peptide_aa_comp, `Amino acid` == ith_aa, Dataset == dataset2)[["Frequency"]])*100, 2)) %>% 
      setNames(c("Amino acid", "pvalue", paste0(dataset1, " mean frequency [%]"), paste0(dataset2, " mean frequency [%]"),
                 paste0(dataset1, " median frequency [%]"), paste0(dataset2, " median frequency [%]")))
  }) %>% bind_rows() 
  
  res <- mutate(test_res,
         "Adjusted p-value" = p.adjust(pvalue, "BH")) %>% 
    select(-pvalue)
  
  write.csv(res, paste0(data_path, "aa_comp_table_", dataset1, "_", dataset2, ".csv"), row.names = FALSE)
  xtable(res, caption = "", label = "", align = "ccccccc",
         digits = rep(2, 7),
         display = c("f", "s", rep("f", 4), "E")) %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "aa_comp_table_", dataset1, "_", dataset2, ".txt"))
}


get_benchmark_results_table <- function(benchmark_res, software) {
  if(software == "anticp") {
    benchmark_res <- arrange(benchmark_res, desc(Software))
  }
  write.csv(benchmark_res, paste0(data_path, "benchmark_results_", software, ".csv"))
  benchmark_res %>% 
    xtable(., caption = "", label = "", align = "cccccccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "benchmark_results_", software, ".txt"))
}


get_ngram_list <- function(imp_ngrams_dat) {
  intersect(intersect(imp_ngrams_dat[["acp_amp"]], 
                      imp_ngrams_dat[["acp_neg"]]),
            imp_ngrams_dat[["amp_neg"]]) %>% 
    decode_ngrams() %>% 
    paste0(., collapse = ", ") %>% 
    writeLines(., paste0(data_path, "ngram_list.txt"))
}
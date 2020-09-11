get_cv_pred_table <- function(cv_mer_performance, cv_peptide_performance) {
  tab <- mutate(cv_mer_performance,
                Layer = "Mer") %>% 
    bind_rows(mutate(cv_peptide_performance,
                     Layer = "Peptide")) %>% 
    pivot_longer(Accuracy:Kappa, names_to = "Measure", values_to = "Value") %>% 
    group_by(Layer, Measure) %>%
    summarise(mean = mean(Value),
              sd = sd(Value)) %>% 
    mutate(value = paste0(round(mean, 2), " (+/-", round(sd, 3), ")")) %>% 
    select(-c(mean, sd)) %>%
    pivot_wider(names_from = Measure, values_from = value)
  write.csv(tab, paste0(data_path, "cv_results.csv"), row.names = FALSE)
    xtable(tab, caption = "", label = "", align = "ccccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
      writeLines(paste0(data_path, "cv_results.txt"))
}

get_datasets_table <- function(acp_train, acp_test, amp_train, amp_test, neg_train, neg_test) {
  df <- data.frame(Dataset = c("Train-test", "Benchmark"),
             ACP = c(length(acp_train), length(acp_test)),
             AMP = c(length(amp_train), length(amp_test)),
             neg = c(length(neg_train), length(neg_test)))
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
    mutate_at(vars(acp, amp, neg), ~round(as.numeric(.), 3))
  colnames(res) <- c("Peptide", "ACP", "AMP", "Negative")
  write.csv(res, paste0(data_path, "mito_ACP_preds.csv"), row.names = FALSE)
  xtable(res, caption = "Fig:mito_ACP_preds", label = "", align = "ccccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "mito_ACP_preds.txt"))
}


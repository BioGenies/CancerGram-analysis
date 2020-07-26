# Select AMPs for negative dataset ensuring the same length distribution as in positive dataset
filter_amps <- function(amp_full_dataset, acp_dataset) {
  amp <- amp_full_dataset[which(!(amp_full_dataset %in% acp_dataset) & lengths(amp_full_dataset) <= 50)] 
  acp_len_groups <- cut(lengths(acp_dataset), 
                        breaks = as.numeric(quantile(lengths(acp_dataset), probs = seq(0, 1, 0.2))),
                        include.lowest = TRUE) %>% 
    table()
  amp_len_groups <- cut(lengths(amp),
                        breaks = c(5, 13, 17, 22, 28, 50),
                        include.lowest = TRUE)
  selected <- sapply(names(acp_len_groups), USE.NAMES = FALSE, function(ith_len_group) {
    n <- acp_len_groups[ith_len_group]
    len_group <- which(amp_len_groups == ith_len_group)
    sample(len_group, n)
  }) %>% unlist()
  amp[selected]
}


get_mers_mc <- function(pos, pos_id, amp, amp_id, neg, neg_id) {
  seq_groups <- lapply(names(pos_id), function(i)
    c(pos[pos_id[[i]][["traintest"]]], amp[amp_id[[i]][["traintest"]]], neg[neg_id[[i]][["traintest"]]])) %>% 
    setNames(names(pos_id))
  
  lapply(names(seq_groups), function(ith_group_id) {
    ith_group <- seq_groups[[ith_group_id]]
    
    folded <- cvFolds(length(ith_group), K = 5)
    fold_df <- data.frame(source_peptide = names(ith_group)[folded[["subsets"]]], 
                          fold = folded[["which"]],
                          stringsAsFactors = FALSE)
    
    ith_group %>% 
      list2matrix() %>% 
      create_mer_df %>% 
      mutate(group = ith_group_id,
             target = case_when(grepl("CUTTED", source_peptide) ~ "neg",
                                grepl("dbAMP", source_peptide) ~ "amp",
                                grepl("CancerPPD|AP|DRAMP", source_peptide) ~ "acp")) %>% 
      inner_join(fold_df, by = c("source_peptide" = "source_peptide"))
  }) %>% 
    do.call(rbind, .)
}  


write_benchmark_mc <- function(pos, pos_id, amp, amp_id, neg, neg_id) {
  seq_list <- c(pos[unlist(lapply(pos_id, function(ith_len_group) ith_len_group[["benchmark"]]))],
                amp[unlist(lapply(amp_id, function(ith_len_group) ith_len_group[["benchmark"]]))],
                neg[unlist(lapply(neg_id, function(ith_len_group) ith_len_group[["benchmark"]]))]) 
  write_fasta(seq_list, file = "results/benchmark_mc.fasta") 
  print(paste0("Number of sequences in the benchmark dataset: ", length(seq_list))) 
}


# 
# get_mc_mer_df <- function(datasets, names) {
#   lapply(1L:length(names), function(i) { 
#     mer_df_from_list(datasets[i]) %>% 
#       mutate(target = names[i]) 
#   }) %>% bind_rows()
# }


get_imp_ngrams_mc <- function(ngrams, mer_df, cutoff = 0.05) {
  combns <- combn(unique(mer_df[["target"]]), 2, simplify = FALSE)
  lapply(combns, function(ith_cmbn) {
    tar <- filter(mer_df, target %in% ith_cmbn) %>% 
      mutate(target = ifelse(target == ith_cmbn[1], 1, 0))
    features <- ngrams[mer_df[["target"]] %in% ith_cmbn,]
    test_bis <- test_features(tar[["target"]], features)
    res <- cut(test_bis, breaks = c(0, cutoff, 1))[[1]]
  }) %>% setNames(sapply(combns, function(ith_combn) paste(ith_combn, collapse = "_")))
}



train_mc_model_mers <- function(mer_df, binary_ngrams, imp_bigrams) {
  ranger_train_data <- data.frame(as.matrix(binary_ngrams[, imp_bigrams]),
                                  tar = as.factor(mer_df[["target"]]))
  model_full_alphabet <- ranger(dependent.variable.name = "tar", data = ranger_train_data, 
                                write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                                verbose = FALSE, classification = TRUE)
  model_full_alphabet
}


calculate_statistics_single <- function(mer_preds, group) {
  if(!("target" %in% colnames(mer_preds))) {
    mer_preds <- mutate(mer_preds, target = case_when(grepl("CUTTED", source_peptide) ~ "neg",
                                                      grepl("dbAMP", source_peptide) ~ "amp",
                                                      grepl("CancerPPD|AP|DRAMP", source_peptide) ~ "acp"))
  }
  mer_preds[c("source_peptide", "target", group)] %>% 
    setNames(c("source_peptide", "target", "pred")) %>% 
    calculate_statistics() %>% {
      df <- .
      nondescriptive_names <- setdiff(colnames(df), c("source_peptide", "target"))
      colnames(df)[colnames(df) %in% nondescriptive_names] <- 
        paste0(group, "_", colnames(df)[colnames(df) %in% nondescriptive_names])
      df
    }
}


train_mc_model_peptides <- function(mer_statistics) {
  train_dat <- mer_statistics %>% 
    select(-source_peptide)
  peptide_model <- ranger(dependent.variable.name = "target", data = train_dat, 
                          write.forest = TRUE, probability = TRUE, num.trees = 500, 
                          verbose = FALSE, classification = TRUE)
  peptide_model
}

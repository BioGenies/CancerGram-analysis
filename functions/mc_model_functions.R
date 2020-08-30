# Select AMPs for negative dataset ensuring the same length distribution as in positive dataset
filter_amps <- function(amp_full_dataset, acp_dataset) {
  amp <- amp_full_dataset[which(!(amp_full_dataset %in% acp_dataset) & lengths(amp_full_dataset) <= 50)] 
  acp_len_groups <- cut(lengths(acp_dataset), 
                        breaks = as.numeric(quantile(lengths(acp_dataset), probs = seq(0, 1, 0.2))),
                        include.lowest = TRUE) %>% 
    table()
  amp_len_groups <- cut(lengths(amp),
                        breaks = as.numeric(quantile(lengths(acp_dataset), probs = seq(0, 1, 0.2))),
                        include.lowest = TRUE)
  selected <- sapply(names(acp_len_groups), USE.NAMES = FALSE, function(ith_len_group) {
    n <- acp_len_groups[ith_len_group]
    len_group <- which(amp_len_groups == ith_len_group)
    sample(len_group, n)
  }) %>% unlist()
  amp[selected]
}


get_mers_mc <- function(acp, acp_id, amp, amp_id, neg, neg_id) {
  seq_groups <- lapply(names(acp_id), function(i)
    c(acp[acp_id[[i]][["traintest"]]], amp[amp_id[[i]][["traintest"]]], neg[neg_id[[i]][["traintest"]]])) %>% 
    setNames(names(acp_id))
  
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
                                grepl("CancerPPD|AP|DRAMP", source_peptide) ~ "acp",
                                grepl("pos_train_main", source_peptide) ~ "acp",
                                grepl("neg_train_main", source_peptide) ~ "amp",
                                grepl("neg_train_alternate", source_peptide) ~ "neg")) %>% 
      inner_join(fold_df, by = c("source_peptide" = "source_peptide"))
  }) %>% 
    do.call(rbind, .)
}  


write_benchmark_mc <- function(acp, acp_id, amp, amp_id, neg, neg_id) {
  seq_list <- c(acp[unlist(lapply(acp_id, function(ith_len_group) ith_len_group[["benchmark"]]))],
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
  
  if ("fold" %in% colnames(mer_preds)) {
    mer_preds[c("source_peptide", "fold", group)] %>% 
      setNames(c("source_peptide", "fold", "pred")) %>% 
      calculate_statistics() %>% {
        df <- .
        nondescriptive_names <- setdiff(colnames(df), c("source_peptide", "target", "fold"))
        colnames(df)[colnames(df) %in% nondescriptive_names] <- 
          paste0(group, "_", colnames(df)[colnames(df) %in% nondescriptive_names])
        df
      } 
  } else {
    mer_preds[c("source_peptide", group)] %>% 
      setNames(c("source_peptide", "pred")) %>% 
      calculate_statistics() %>% {
        df <- .
        nondescriptive_names <- setdiff(colnames(df), c("source_peptide", "target"))
        colnames(df)[colnames(df) %in% nondescriptive_names] <- 
          paste0(group, "_", colnames(df)[colnames(df) %in% nondescriptive_names])
        df
      }
  }
}


calculate_statistics_mc <- function(mer_preds, groups) {
  res <- lapply(groups, function(i) {
    calculate_statistics_single(mer_preds, i)
  }) %>% do.call(cbind, .) 
  stat <- res[,!duplicated(colnames(res))] %>% 
    select(-c("amp_n_peptide", "neg_n_peptide"))
  stat
} 


train_mc_model_peptides <- function(mer_statistics) {
  if("fold" %in% colnames(mer_statistics)) {
    train_dat <- mer_statistics %>% 
      select(-c(source_peptide, fold))
  } else {
    train_dat <- mer_statistics %>% 
      select(-source_peptide)
  }
  peptide_model <- ranger(dependent.variable.name = "target", data = train_dat, 
                          write.forest = TRUE, probability = TRUE, num.trees = 500, 
                          verbose = FALSE, classification = TRUE)
  peptide_model
}


do_cv_mc <- function(mer_df, binary_ngrams, cutoff) {
  lapply(unique(mer_df[["fold"]]), function(ith_fold) {
    print(paste0(ith_fold))
    train_dat <- filter(mer_df, fold != ith_fold)
    test_dat <- filter(mer_df, fold == ith_fold)
    
    ngrams_to_test <- binary_ngrams[mer_df[["fold"]] != ith_fold, ]
    
    imp_bigrams <- unique(unlist(unname(get_imp_ngrams_mc(ngrams_to_test, train_dat, cutoff))))
    
    ranger_train_data <- data.frame(as.matrix(binary_ngrams[mer_df[["fold"]] != ith_fold, imp_bigrams]),
                                    tar = as.factor(train_dat[["target"]]))
    model_cv <- ranger(dependent.variable.name = "tar", data =  ranger_train_data, 
                       write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                       verbose = FALSE)
    
    preds <- cbind(test_dat,
                   predict(model_cv, 
                           data.frame(as.matrix(binary_ngrams[mer_df[["fold"]] == ith_fold, imp_bigrams])))[["predictions"]])
    
    preds[, c("source_peptide", "mer_id", 
              "fold", "target", "acp", "amp", "neg")] 
  }) %>% bind_rows()
}

do_cv_peptides_mc <- function(cv_mer_res) {
  lapply(unique(cv_mer_res[["fold"]]), function(ith_fold) {
    print(paste0(ith_fold))
    cv_stats <- cv_mer_res %>% 
      calculate_statistics_mc(c("acp", "amp", "neg")) %>% 
      mutate(target = factor(case_when(grepl("pos_train_main|Cancer|AP|DRAMP", source_peptide) ~ "acp",
                                       grepl("neg_train_main|dbAMP", source_peptide) ~ "amp",
                                       grepl("neg_train_alternate|CUTTED", source_peptide) ~ "neg")))
    
    train_dat <- filter(cv_stats, fold != ith_fold)
    test_dat <- filter(cv_stats, fold == ith_fold)
    
    model_cv <- train_mc_model_peptides(train_dat)
    
    preds <- cbind(test_dat[, c("source_peptide", "target", "fold")],
                   predict(model_cv, test_dat)[["predictions"]])
    
    preds[, c("source_peptide", "fold", "target", "acp", "amp", "neg")] 
    
  }) %>% bind_rows()
}




calc_cv_performance <- function(preds) {
  lapply(unique(preds[["fold"]]), function(ith_fold){
    acp_amp <- filter(preds, grepl("train_main", source_peptide) & fold == ith_fold) %>% 
      mutate(decision = factor(ifelse((acp > amp & acp > neg), "TRUE", "FALSE")),
             target = factor(ifelse(target == "acp", "TRUE", "FALSE")))
    acp_neg <- filter(preds, grepl("pos_train_main|neg_train_alt", source_peptide) & fold == ith_fold) %>% 
      mutate(decision = factor(ifelse((acp > amp & acp > neg), "TRUE", "FALSE")),
             target = factor(ifelse(target == "acp", "TRUE", "FALSE")))
    datasets <- list("ACP/AMP" = acp_amp, "ACP/neg" = acp_neg)
    lapply(names(datasets), function(ith_set) {
      predictions <- datasets[[ith_set]]
      data.frame(
        dataset = ith_set,
        fold = ith_fold,
        AUC = mlr3measures::auc(predictions[["target"]], predictions[["acp"]], "TRUE"),
        MCC = mlr3measures::mcc(predictions[["target"]], predictions[["decision"]], "TRUE"),
        Precision = mlr3measures::precision(predictions[["target"]], predictions[["decision"]], "TRUE"),
        Sensitivity = mlr3measures::sensitivity(predictions[["target"]], predictions[["decision"]], "TRUE"),
        Specificity = mlr3measures::specificity(predictions[["target"]], predictions[["decision"]], "TRUE"),
        Accuracy = mlr3measures::acc(predictions[["target"]], predictions[["decision"]]),
        stringsAsFactors = FALSE)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

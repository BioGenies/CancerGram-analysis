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

train_model_mers <- function(train_mer_df, train_binary_ngrams, imp_bigrams) {
  ranger_train_data <- data.frame(as.matrix(train_binary_ngrams[, imp_bigrams]),
                                  tar = as.factor(train_mer_df[["target"]]))
  model_mer <- ranger(dependent.variable.name = "tar", data = ranger_train_data, 
                      write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                      verbose = FALSE)
  model_mer
}

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)),
                   na.last = FALSE)
  x[x_order]
}



train_and_test_anticp <- function(pos_train, neg_train, pos_test, neg_test) {
  
  mer_df <- lapply(list(pos_train, neg_train, pos_test, neg_test), function(ith_set) {
    mer_df_from_list(ith_set) %>% 
      mutate(target = ifelse(grepl("pos", source_peptide), "TRUE", "FALSE"),
             dataset = ifelse(grepl("train", source_peptide), "train", "test")) 
  }) %>% bind_rows()
  
  ngrams <- count_and_gather_ngrams(mer_df,
                                    c(1, rep(2, 4), rep(3, 4)),
                                    list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1)))
  
  train_dat <- filter(mer_df, dataset == "train")
  test_dat <- filter(mer_df, dataset == "test")
  
  test_bis <- test_features(as.logical(train_dat[["target"]]),
                            ngrams[mer_df[["dataset"]] == "train", ])
  
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  ranger_train_data <- data.frame(as.matrix(ngrams[mer_df[["dataset"]] == "train", imp_bigrams]),
                                  tar = as.factor(train_dat[["target"]]))
  model <- ranger(dependent.variable.name = "tar", data =  ranger_train_data, 
                  write.forest = TRUE, probability = TRUE, num.trees = 500, 
                  verbose = FALSE, seed = 2938)
  
  preds <- mutate(test_dat,
                  pred = predict(model, 
                                 data.frame(as.matrix(ngrams[mer_df[["dataset"]] == "test", imp_bigrams])))[["predictions"]][, "TRUE"])
}





get_metrics <- function(preds, cutoff, dataset_name) {
  res <- preds %>% 
    select(c(source_peptide, target, pred)) %>% 
    group_by(source_peptide) %>% 
    summarise(decision = ifelse(any(pred > cutoff), TRUE, FALSE)) %>% 
    inner_join(preds[, c("source_peptide", "target")]) %>% 
    unique() %>% 
    mutate(target = factor(target),
           decision = factor(decision))
  
  data.frame(MCC = mlr3measures::mcc(res[["target"]], res[["decision"]], "TRUE"),
             Precision = mlr3measures::precision(res[["target"]], res[["decision"]], "TRUE"),
             Sensitivity = mlr3measures::sensitivity(res[["target"]], res[["decision"]], "TRUE"),
             Specificity = mlr3measures::specificity(res[["target"]], res[["decision"]], "TRUE"),
             Accuracy = mlr3measures::acc(res[["target"]], res[["decision"]]),
             stringsAsFactors = FALSE) %>% 
    mutate(dataset = dataset_name,
           cutoff = cutoff)
}




get_validation_dataset <- function() {
  
  benchmark_our <- read_fasta("./results/benchmark_mc.fasta")
  benchmark_anticp <- c(neg_test_main, neg_test_alt, pos_test_main)
  
  traintest_acp <- sapply(names(cdhit_acp_data_ids), function(x) cdhit_acp_data[cdhit_acp_data_ids[[x]][["benchmark"]]])
  traintest_amp <- sapply(names(amp_filtered_data_ids), function(x) amp_filtered_data[amp_filtered_data_ids[[x]][["benchmark"]]])
  traintest_neg <- sapply(names(neg_data_ids), function(x) neg_data[neg_data_ids[[x]][["benchmark"]]]) 
  
  traintest_our <- unlist(c(traintest_acp, traintest_amp, traintest_neg), recursive = FALSE) 
  traintest_anticp <- c(neg_train_main, neg_train_alt, pos_train_main)
  
  c(benchmark_anticp[which(!(benchmark_anticp %in% traintest_our))],
    benchmark_our[which(!(benchmark_our %in% traintest_anticp))])
}



calc_performance <- function(preds) {
  acp_amp <- filter(preds, grepl("test_main|train_main", source_peptide)) %>% 
    mutate(decision = factor(ifelse(acp >= 0.5, TRUE, FALSE)),
           target = factor(ifelse(target == "acp", TRUE, FALSE)))
  acp_neg <- filter(preds, grepl("test_alt|train_alt", source_peptide)) %>% 
    mutate(decision = factor(ifelse(acp >= 0.5, TRUE, FALSE)),
           target = factor(ifelse(target == "acp", TRUE, FALSE)))
  datasets <- list("ACP/AMP" = acp_amp, "ACP/neg" = acp_neg)
  lapply(names(datasets), function(ith_set) {
    predictions <- datasets[[ith_set]]
    data.frame(
      dataset = ith_set,
      AUC = mlr3measures::auc(predictions[["target"]], predictions[["acp"]], "TRUE"),
      MCC = mlr3measures::mcc(predictions[["target"]], predictions[["decision"]], "TRUE"),
      Precision = mlr3measures::precision(predictions[["target"]], predictions[["decision"]], "TRUE"),
      Sensitivity = mlr3measures::sensitivity(predictions[["target"]], predictions[["decision"]], "TRUE"),
      Specificity = mlr3measures::specificity(predictions[["target"]], predictions[["decision"]], "TRUE"),
      Accuracy = mlr3measures::acc(predictions[["target"]], predictions[["decision"]]),
      stringsAsFactors = FALSE)
  }) %>% bind_rows()
}


write_benchmark_mc <- function(acp, acp_id, amp, amp_id, neg, neg_id) {
  seq_list <- c(acp[unlist(lapply(acp_id, function(ith_len_group) ith_len_group[["benchmark"]]))],
                amp[unlist(lapply(amp_id, function(ith_len_group) ith_len_group[["benchmark"]]))],
                neg[unlist(lapply(neg_id, function(ith_len_group) ith_len_group[["benchmark"]]))]) 
  write_fasta(seq_list, file = "results/benchmark_mc.fasta") 
  print(paste0("Number of sequences in the benchmark dataset: ", length(seq_list))) 
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


#' Get dataframe of mers
#' 
#' This function creates a dataframe of 10-mers with columns indicating 
#' their source peptide, mer ID, group of source peptide, fold and target.
#' @param pos positive dataset
#' @param pos_id IDs of the positive dataset
#' @param neg negative dataset
#' @param neg_id IDs of the negative dataset
#' @return a dataframe of 10-mers with their source peptides, mer IDs, group,
#' target and fold
get_mers <- function(pos, pos_id, neg, neg_id) {
  seq_groups <- lapply(names(pos_id), function(i)
    c(pos[pos_id[[i]][["traintest"]]], neg[neg_id[[i]][["traintest"]]])) %>% 
    setNames(names(pos_id))
  #lapply(pos_id, function(i) i[["traintest"]])
  #lapply(neg_id, function(i) i[["traintest"]])
  
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
             target = ifelse(grepl("CUTTED", source_peptide), FALSE, TRUE)) %>% 
      inner_join(fold_df, by = c("source_peptide" = "source_peptide"))
  }) %>% 
    do.call(rbind, .)
}  


calc_imp_bigrams <- function(train_mer_df, train_binary_ngrams, cutoff = 0.05) {
  test_bis <- test_features(train_mer_df[["target"]], train_binary_ngrams)
  imp_bigrams <- cut(test_bis, breaks = c(0, cutoff, 1))[[1]]
  imp_bigrams
}



train_model_peptides <- function(mer_statistics) {
  train_dat <- mer_statistics %>% 
    select(c("target", "fraction_true", "pred_mean", "pred_median",
             "n_peptide", "n_pos", "pred_min", "pred_max", "longest_pos",
             "n_pos_5", "frac_0_0.2", "frac_0.2_0.4", "frac_0.4_0.6",
             "frac_0.6_0.8", "frac_0.8_1"))
  model_cv <- ranger(dependent.variable.name = "target", data = train_dat, 
                     write.forest = TRUE, probability = TRUE, num.trees = 500, 
                     verbose = FALSE, classification = TRUE)
  model_cv
}



do_cv <- function(mer_df, binary_ngrams, cutoff = 0.05) {
  possible_groups <- sort_group(unique(mer_df[["group"]]))
  group_combs <- list(possible_groups[1:2], possible_groups[1:3], possible_groups[1:4], c(possible_groups))
  lapply(group_combs, function(ith_group) {
    lapply(unique(mer_df[["fold"]]), function(ith_fold) {
      print(paste0(ith_group, "|", ith_fold))
      train_dat <- filter(mer_df, group %in% ith_group, fold != ith_fold)
      test_dat <- filter(mer_df, fold == ith_fold)
      
      test_bis <- test_features(train_dat[["target"]],
                                binary_ngrams[mer_df[["group"]] %in% ith_group & 
                                                mer_df[["fold"]] != ith_fold, ])
      
      imp_bigrams <- cut(test_bis, breaks = c(0, cutoff, 1))[[1]]
      
      ranger_train_data <- data.frame(as.matrix(binary_ngrams[mer_df[["group"]] %in% ith_group & 
                                                                mer_df[["fold"]] != ith_fold, imp_bigrams]),
                                      tar = as.factor(train_dat[["target"]]))
      model_cv <- ranger(dependent.variable.name = "tar", data =  ranger_train_data, 
                         write.forest = TRUE, probability = TRUE, num.trees = 500, 
                         verbose = FALSE)
      
      preds <- mutate(test_dat,
                      pred = predict(model_cv, 
                                     data.frame(as.matrix(binary_ngrams[mer_df[["fold"]] == ith_fold, imp_bigrams])))[["predictions"]][, "TRUE"],
                      comb = paste(ith_group, collapse = ","))
      
      # single mer predictions
      #HMeasure(preds[["target"]], preds[["pred"]])[["metrics"]]
      
      preds[, c("source_peptide", "mer_id", "group", 
                "fold", "target", "pred", "comb")] 
    }) %>% bind_rows()
  }) %>% bind_rows()
}


do_cv_degenerate <- function(mer_df, binary_ngrams, elements_groups, cutoff, mc = TRUE) {
  pblapply(elements_groups, function(ith_alphabet){
    deg_binary_ngrams <- degenerate_ngrams(binary_ngrams, string2list(ith_alphabet), binarize = TRUE)
    (if (mc == TRUE) {
      do_cv_mc(mer_df, deg_binary_ngrams, cutoff)
    } else {
      do_cv(mer_df, deg_binary_ngrams, cutoff)
    })%>% 
      mutate(alphabet = ith_alphabet) %>% 
      write.csv(file = paste0(data_path, "results/", ith_alphabet, ".csv"), 
                row.names = FALSE)
  })
}

sort_group <- function(x) {
  splitted_x <- sapply(strsplit(x, split = ","), function(i) i[1])
  x_order <- order(as.numeric(gsub(pattern = "[^0-9]", 
                                   replacement = "", x = splitted_x)))
  x[x_order]
}

string2list <- function(x) {
  pasted_group <- strsplit(x, "_", fixed = TRUE)[[1]] %>% 
    toupper()
  res <- strsplit(pasted_group, "")
  names(res) <- 1L:length(res)
  res
}

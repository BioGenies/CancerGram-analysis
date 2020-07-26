calc_imp_bigrams <- function(train_mer_df, train_binary_ngrams, cutoff = 0.05) {
  test_bis <- test_features(train_mer_df[["target"]], train_binary_ngrams)
  imp_bigrams <- cut(test_bis, breaks = c(0, cutoff, 1))[[1]]
  imp_bigrams
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

count_imp_ngrams <- function(mer_df, imp_ngrams) {
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_specified(imp_ngrams) %>% 
    binarize 
}


get_single_seq_mers <- function(seq) {
  seq2ngrams(seq, 10, a()[-1]) %>% 
    decode_ngrams() %>% 
    unname() %>% 
    strsplit(split = "") %>% 
    do.call(rbind, .) 
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

# Processing AntiCP sequence files
process_sequences <- function(seq_file) {
  seqs <- readLines(paste0("./data/", seq_file))
  named_seqs <- seqs %>% 
    strsplit("") %>% 
    setNames(paste0(strsplit(seq_file, ".", fixed = TRUE)[[1]][1], "_", 1:length(seqs)))
  named_seqs[lengths(named_seqs) >= 5]
}
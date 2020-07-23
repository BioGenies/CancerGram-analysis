library(dplyr)
library(biogram)
library(seqR)
library(ranger)

source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")

process_sequences <- function(seq_file) {
  seqs <- readLines(paste0("./data/", seq_file))
  named_seqs <- seqs %>% 
    strsplit("") %>% 
    setNames(paste0(strsplit(seq_file, ".", fixed = TRUE)[[1]][1], "_", 1:length(seqs)))
  named_seqs[lengths(named_seqs) >= 5]
}


train_and_test <- function(pos_train, neg_train, pos_test, neg_test) {
  mer_df <- lapply(list(pos_train, neg_train, pos_test, neg_test), function(ith_set) {
    ith_set %>% 
      list2matrix() %>% 
      create_mer_df() %>% 
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


# Read main dataset sequences (neg = random peptides)
pos_train_main <- process_sequences("pos_train_main.txt") # 3 sequences < 5 aa
pos_test_main <- process_sequences("pos_test_main.txt") # 1 sequence < 5 aa
neg_train_main <- process_sequences("neg_train_main.txt") 
neg_test_main <- process_sequences("neg_test_main.txt") # 2 sequences < 5 aa

# Read alternate dataset sequences (neg = AMPs)
pos_train_alt <- process_sequences("pos_train_alternate.txt") # 4 sequences < 5 aa
pos_test_alt <- process_sequences("pos_test_alternate.txt") 
neg_train_alt <- process_sequences("neg_train_alternate.txt")
neg_test_alt <- process_sequences("neg_test_alternate.txt")


# Get predictions 
main_preds <- train_and_test(pos_train_main, neg_train_main, pos_test_main, neg_test_main)
alt_preds <- train_and_test(pos_train_alt, neg_train_alt, pos_test_alt, neg_test_alt)

# Calculate metrics
datasets <- list(main_preds, alt_preds)
cutoffs <- c(rep(0.5, 2), rep(0.7, 2), rep(0.9, 2))
dataset_names <- c("main", "alt")

result <- mapply(function(set, cutoff, name) {
  get_metrics(set, cutoff, name)
}, set = datasets, cutoff = cutoffs, name = dataset_names)

saveRDS(main_preds, "./results/main_dataset_preds.RDS")
saveRDS(alt_preds, "./results/alt_dataset_preds.RDS")

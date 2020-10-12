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

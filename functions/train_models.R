train_mc_model_mers <- function(mer_df, binary_ngrams, imp_bigrams) {
  ranger_train_data <- data.frame(as.matrix(binary_ngrams[, imp_bigrams]),
                                  tar = as.factor(mer_df[["target"]]))
  model_full_alphabet <- ranger(dependent.variable.name = "tar", data = ranger_train_data, 
                                write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                                verbose = FALSE, classification = TRUE)
  model_full_alphabet
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

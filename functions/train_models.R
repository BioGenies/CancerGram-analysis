#' Train mer model
#' 
#' Train a multiclass model on mer level.
#' @param mer_df data frame of mers
#' @param binary_ngrams binarized n-gram occurences
#' @param imp_bigrams \code{character} vector of informative n-grams
#' @return random forest model of class \code{ranger}
train_mc_model_mers <- function(mer_df, binary_ngrams, imp_bigrams) {
  ranger_train_data <- data.frame(as.matrix(binary_ngrams[, imp_bigrams]),
                                  tar = as.factor(mer_df[["target"]]))
  model_full_alphabet <- ranger(dependent.variable.name = "tar", data = ranger_train_data, 
                                write.forest = TRUE, probability = TRUE, num.trees = 2000, 
                                verbose = FALSE, classification = TRUE)
  model_full_alphabet
}


#' Train peptide model
#' 
#' Train a multiclass model on peptide level.
#' @param mer_statistics data frame of statistics calculated for mer 
#' predictions using \link{\code{calculate_statistics_mc}}
#' @return random forest model of class \code{ranger}
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

#' Get important ngrams from multiclass comparisons
#' 
#' This function performs feature selection using QuiPT for all
#' binary combinations of target classes.
#' @param ngrams binarized n-gram counts
#' @param mer_df data frame of mers
#' @param cutoff p-value treshold for QuiPT
#' @return list of length equal to all possible combinations
#' of target classes. Each element of a list is a \code{character}
#' vector of informative n-grams for a comparison of two classes.
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

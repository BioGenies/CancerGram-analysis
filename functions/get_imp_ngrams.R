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

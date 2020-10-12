count_ngrams <- function(mer_df, k, gaps) {
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_kmers(sequences = .,
               k = k,
               kmer_gaps = gaps,
               alphabet = toupper(colnames(aaprop))) %>% 
    binarize
}


count_and_gather_ngrams <- function(mer_df, k_list, gap_list) {
  res <- mapply(function(k, gap) {
    count_ngrams(mer_df, k, gap)
  }, k = k_list, gap = gap_list, SIMPLIFY = FALSE) %>% 
    do.call(cbind, .)
  colnames(res)[which(nchar(colnames(res)) == 1)] <- paste0(colnames(res)[which(nchar(colnames(res)) == 1)], "_0")
  res
}



count_imp_ngrams <- function(mer_df, imp_ngrams) {
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_specified(imp_ngrams) %>% 
    binarize 
}

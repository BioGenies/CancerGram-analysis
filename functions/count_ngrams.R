#' Count n-grams
#' 
#' Counts occurences of n-grams in data frame of mers and then
#' binarizes them.
#' @param mer_df data frame of mers
#' @param k length of n-grams to count
#' @param gaps \code{numeric} vector of gaps sizes
#' @return \code{simple triplet matrix} with n-grams as columns
#' and their binarized occurences in mers as rows
count_ngrams <- function(mer_df, k, gaps) {
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_kmers(sequences = .,
               k = k,
               kmer_gaps = gaps,
               alphabet = toupper(colnames(aaprop))) %>% 
    binarize
}

#' Count and gather n-grams
#' 
#' Counts occurences of n-grams in a data frame of mers, binarizes them 
#' and combines results for n-grams of different sizes.
#' @param mer_df data frame of mers
#' @param k_list \code{list} of lengths of n-grams to count
#' @param gap_list \code{list} of \code{numeric} vectors of gaps sizes
#' @return \code{simple triplet matrix} with n-grams as columns
#' and their binarized occurences in mers as rows
count_and_gather_ngrams <- function(mer_df, k_list, gap_list) {
  res <- mapply(function(k, gap) {
    count_ngrams(mer_df, k, gap)
  }, k = k_list, gap = gap_list, SIMPLIFY = FALSE) %>% 
    do.call(cbind, .)
  colnames(res)[which(nchar(colnames(res)) == 1)] <- paste0(colnames(res)[which(nchar(colnames(res)) == 1)], "_0")
  res
}


#' Count important n-grams
#' 
#' Counts occurences of specified informative n-grams in a data frame 
#' of mers and binarizes them.
#' @param mer_df data frame of mers
#' @param imp_ngrams \code{character} vector of informative n-grams to count
#' @return \code{matrix} of binarized n-gram occurences in mers
count_imp_ngrams <- function(mer_df, imp_ngrams) {
  mer_df[, grep("^X", colnames(mer_df))] %>% 
    as.matrix() %>% 
    count_specified(imp_ngrams) %>% 
    binarize 
}

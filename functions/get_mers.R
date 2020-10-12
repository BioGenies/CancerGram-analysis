#' Create mer dataframe
#' 
#' Creates dataframe of 5-mers
#' @param seq matrix of sequences created by \code{\link{list2matrix}}
#' @return dataframe of mers with their source peptide and mer ID
create_mer_df <- function(seq) 
  do.call(rbind, lapply(1L:nrow(seq), function(i) {
    seq2ngrams(seq[i, ][!is.na(seq[i, ])], 5, a()[-1]) %>% 
      decode_ngrams() %>% 
      unname() %>% 
      strsplit(split = "") %>% 
      do.call(rbind, .) %>% 
      data.frame(stringsAsFactors = FALSE) %>% 
      mutate(source_peptide = rownames(seq)[i],
             mer_id = paste0(source_peptide, "m", 1L:nrow(.)))
  }))


#' Get mer data frame from sequence list
#' 
#' Creates data frame of mers from a list of sequences.
#' @param seq_list list of sequences
#' @return data frame of mers with their source peptide and mer ID
mer_df_from_list <- function(seq_list) {
  seq_list %>% 
    list2matrix() %>% 
    create_mer_df() 
}


#' Get mer data frame with length groups
#' 
#' Creates data frame of mers from a list of list of sequences, 
#' divides the sequences into equally distributed length groups 
#' and based on those groups performs classification into 5 folds.
#' @param list_of_seq_list list of sequence lists
#' @return data frame of mers with their source peptide, mer ID,
#' length group and fold.
mer_df_from_list_len_group <- function(list_of_seq_list) {
  lapply(list_of_seq_list, function(ith_seq_list) {
    
    lens <- data.frame(source_peptide = names(ith_seq_list),
                       len_group = cut(lengths(ith_seq_list),
                                       breaks = as.numeric(quantile(lengths(ith_seq_list), probs = seq(0, 1, 0.2))),
                                       include.lowest = TRUE)) 
    
    lapply(unique(lens[["len_group"]]), function(ith_group_id) {
      ith_group <- filter(lens, len_group == ith_group_id)[["source_peptide"]]
      
      folded <- cvFolds(length(ith_group), K = 5)
      fold_df <- data.frame(source_peptide = ith_group[folded[["subsets"]]], 
                            fold = folded[["which"]],
                            stringsAsFactors = FALSE)
      
      mer_df <- left_join(mer_df_from_list(ith_seq_list),
                          lens, 
                          by = "source_peptide") %>% 
        inner_join(fold_df, by = c("source_peptide" = "source_peptide"))
    }) %>% 
      do.call(rbind, .)
  }) %>% bind_rows()
}


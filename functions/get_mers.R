#' Create mer dataframe
#' 
#' Creates dataframe of 10-mers
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



mer_df_from_list <- function(seq_list) {
  seq_list %>% 
    list2matrix() %>% 
    create_mer_df() 
}



mer_df_from_list_len_group <- function(seq_list) {
  lens <- data.frame(source_peptide = names(seq_list),
                     len_group = cut(lengths(seq_list),
                                     breaks = c(5, 13, 17, 22, 28, 50),
                                     include.lowest = TRUE)) 
  
  lapply(unique(lens[["len_group"]]), function(ith_group_id) {
    ith_group <- filter(lens, len_group == ith_group_id)[["source_peptide"]]
    
    folded <- cvFolds(length(ith_group), K = 5)
    fold_df <- data.frame(source_peptide = ith_group[folded[["subsets"]]], 
                          fold = folded[["which"]],
                          stringsAsFactors = FALSE)
    
    mer_df <- left_join(mer_df_from_list(seq_list),
                        lens, 
                        by = "source_peptide") %>% 
      inner_join(fold_df, by = c("source_peptide" = "source_peptide"))
  }) %>% 
    do.call(rbind, .)
}


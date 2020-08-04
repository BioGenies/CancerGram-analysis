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

#' Get dataframe of mers
#' 
#' This function creates a dataframe of 10-mers with columns indicating 
#' their source peptide, mer ID, group of source peptide, fold and target.
#' @param pos positive dataset
#' @param pos_id IDs of the positive dataset
#' @param neg negative dataset
#' @param neg_id IDs of the negative dataset
#' @return a dataframe of 10-mers with their source peptides, mer IDs, group,
#' target and fold
get_mers <- function(pos, pos_id, neg, neg_id) {
  seq_groups <- lapply(names(pos_id), function(i)
    c(pos[pos_id[[i]][["traintest"]]], neg[neg_id[[i]][["traintest"]]])) %>% 
    setNames(names(pos_id))
  #lapply(pos_id, function(i) i[["traintest"]])
  #lapply(neg_id, function(i) i[["traintest"]])
  
  lapply(names(seq_groups), function(ith_group_id) {
    ith_group <- seq_groups[[ith_group_id]]
    
    folded <- cvFolds(length(ith_group), K = 5)
    fold_df <- data.frame(source_peptide = names(ith_group)[folded[["subsets"]]], 
                          fold = folded[["which"]],
                          stringsAsFactors = FALSE)
    
    ith_group %>% 
      list2matrix() %>% 
      create_mer_df %>% 
      mutate(group = ith_group_id,
             target = ifelse(grepl("CUTTED", source_peptide), FALSE, TRUE)) %>% 
      inner_join(fold_df, by = c("source_peptide" = "source_peptide"))
  }) %>% 
    do.call(rbind, .)
}  


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

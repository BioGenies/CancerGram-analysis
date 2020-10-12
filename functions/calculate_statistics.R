count_longest <- function(x) {
  splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""),
                         split = "0")[[1]]
  len <- unname(sapply(splitted_x, nchar))
  if (length(len[len > 0]) == 0) {
    0 } else {
      len[len > 0]
    }
}

calculate_statistics <- function(pred_mers) {
  (if("fold" %in% colnames(pred_mers)) {
    group_by(pred_mers, source_peptide, fold)
  } else {  
    group_by(pred_mers, source_peptide)
  }) %>% 
    summarise(fraction_true = mean(pred > 0.5),
              pred_mean = mean(pred),
              pred_median = median(pred),
              n_peptide = length(pred),
              n_pos = sum(pred > 0.5),
              pred_min = min(pred),
              pred_max = max(pred), 
              longest_pos = max(count_longest(pred)),
              n_pos_5 = sum(count_longest(pred) >= 5),
              frac_0_0.2 = sum(pred <= 0.2)/n(),
              frac_0.2_0.4 = sum(pred > 0.2 & pred <= 0.4)/n(),
              frac_0.4_0.6 = sum(pred > 0.4 & pred <= 0.6)/n(),
              frac_0.6_0.8 = sum(pred > 0.6 & pred <= 0.8)/n(),
              frac_0.8_1 = sum(pred > 0.8 & pred <= 1)/n()) %>% 
    ungroup() 
}


calculate_statistics_single <- function(mer_preds, group) {
  
  if ("fold" %in% colnames(mer_preds)) {
    mer_preds[c("source_peptide", "fold", group)] %>% 
      setNames(c("source_peptide", "fold", "pred")) %>% 
      calculate_statistics() %>% {
        df <- .
        nondescriptive_names <- setdiff(colnames(df), c("source_peptide", "target", "fold"))
        colnames(df)[colnames(df) %in% nondescriptive_names] <- 
          paste0(group, "_", colnames(df)[colnames(df) %in% nondescriptive_names])
        df
      } 
  } else {
    mer_preds[c("source_peptide", group)] %>% 
      setNames(c("source_peptide", "pred")) %>% 
      calculate_statistics() %>% {
        df <- .
        nondescriptive_names <- setdiff(colnames(df), c("source_peptide", "target"))
        colnames(df)[colnames(df) %in% nondescriptive_names] <- 
          paste0(group, "_", colnames(df)[colnames(df) %in% nondescriptive_names])
        df
      }
  }
}


calculate_statistics_mc <- function(mer_preds, groups) {
  res <- lapply(groups, function(i) {
    calculate_statistics_single(mer_preds, i)
  }) %>% do.call(cbind, .) 
  stat <- res[,!duplicated(colnames(res))] %>% 
    select(-c("amp_n_peptide", "neg_n_peptide"))
  stat
} 

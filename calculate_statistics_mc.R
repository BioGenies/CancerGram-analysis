library(dplyr)

calculate_statistics <- function(pred_mers) {
  (if("fold" %in% colnames(pred_mers)) {
    group_by(pred_mers, source_peptide, target, fold)
  } else {  
    group_by(pred_mers, source_peptide, target)
  }) %>% 
    summarise(fraction_true = mean(pred > 0.5),
              pred_mean = mean(pred),
              pred_median = median(pred),
              n_peptide = length(pred),
              n_pos = sum(pred > 0.5),
              pred_min = min(pred),
              pred_max = max(pred), 
              #longest_pos = max(count_longest(pred)),
              #n_pos_10 = sum(count_longest(pred) >= 10),
              frac_0_0.2 = sum(pred <= 0.2)/n(),
              frac_0.2_0.4 = sum(pred > 0.2 & pred <= 0.4)/n(),
              frac_0.4_0.6 = sum(pred > 0.4 & pred <= 0.6)/n(),
              frac_0.6_0.8 = sum(pred > 0.6 & pred <= 0.8)/n(),
              frac_0.8_1 = sum(pred > 0.8 & pred <= 1)/n()) %>% 
    ungroup() %>% 
    mutate(target = factor(target))
}


df <- data.frame(source_peptide = "pep1", 
                 target = TRUE,
                 name = c("mer1", "mer2"), 
                 acp = c(0.9, 0.6), 
                 amp = c(0.09, 0.3), 
                 neg = c(0.01, 0.1))

calculate_statistics_single <- function(pred_mer, what) {
  pred_mer[c("source_peptide", "target", what)] %>% 
    setNames(c("source_peptide", "target", "pred")) %>% 
    calculate_statistics() %>% {
      df <- .
      nondescriptive_names <- setdiff(colnames(df), c("source_peptide", 
                                      "target", "fold"))
      colnames(df)[colnames(df) %in% nondescriptive_names] <- 
        paste0(what, "_", colnames(df)[colnames(df) %in% nondescriptive_names])
      df
    }
}

calculate_statistics_single(df, "acp")

count_longest <- function(x) {
  splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""),
                         split = "0")[[1]]
  len <- unname(sapply(splitted_x, nchar))
  if (length(len[len > 0]) == 0) {
    0 } else {
      len[len > 0]
    }
}


test_alphabet_mc <- function(alphabet_file) {
    dat <- read.csv(alphabet_file, stringsAsFactors = FALSE)
    stats <- calculate_statistics_mc(dat, c("acp", "amp", "neg")) %>% 
      mutate(target = factor(case_when(grepl("pos_train_main", source_peptide) ~ "acp",
                                       grepl("neg_train_main", source_peptide) ~ "amp",
                                       grepl("neg_train_alternate", source_peptide) ~ "neg")))
    lapply(1L:5, function(ith_fold) {
      test_dat <- filter(stats, fold == ith_fold) %>% 
        select(-c("source_peptide", "target", "fold"))
      trained_model <- stats %>% 
        filter(fold != ith_fold) %>% 
        select(-fold) %>% 
        train_mc_model_peptides()
      preds <- cbind(filter(stats, fold == ith_fold)[, c("source_peptide", "target", "fold")],
                     predict(trained_model,
                             test_dat)[["predictions"]]) %>% 
        mutate(alphabet = dat[["alphabet"]][1])
    }) %>% bind_rows()
}


test_all_alphabets_mc <- function(data_path, alphabets) {
  lapply(alphabets, function(ith_alphabet) {
    test_alphabet_mc(paste0(data_path, "results/", ith_alphabet, ".csv"))
  }) %>% bind_rows()
}
  

calc_measures_alphabets_mc <- function(alphabets_preds) {
  lapply(unique(alphabets_preds[["alphabet"]]), function(ith_alphabet) {
    lapply(unique(alphabets_preds[["fold"]]), function(ith_fold) {
      dat <- filter(alphabets_preds, fold == ith_fold & alphabet == ith_alphabet) %>% 
        mutate(decision = as.factor(case_when(acp > amp & acp > neg  ~ "acp",
                                              amp > acp & amp > neg ~ "amp",
                                              neg > amp & neg > acp ~ "neg")))
      data.frame(alphabet = ith_alphabet,
                 fold = ith_fold,
                 AU1U = multiclass.AU1U(dat[, c("acp", "amp", "neg")], dat[["target"]]),
                 kappa = KAPPA(dat[["target"]], dat[["decision"]]))
    }) %>% bind_rows()
  }) %>% bind_rows()
}

library(dplyr)
library(drake)
library(biogram)
library(magrittr)
library(seqR)
library(ranger)
library(ggplot2)
library(tidyr)
library(cvTools)
library(measures)

source("functions/mc_model_functions.R")
source("functions/do_cv.R")
source("functions/train_model_peptides.R")
source("functions/benchmark_functions.R")
source("functions/publication_functions.R")
source("functions/get_mers.R")


loadd(c(mers_mc_anticp, ngrams_mc_anticp,
        pos_train_main, pos_test_main, 
        neg_train_main, neg_test_main, neg_train_alt, neg_test_alt), 
      path = "benchmark_cache")


cutoffs <- c(0.001, 0.0001, 1e-5, 1e-8, 1e-10, 1e-15)

test_cutoffs_cv <- function(cutoffs) {
  lapply(cutoffs, function(ith_cutoff) {
    mer_cv <- do_cv_mc(mers_mc_anticp, ngrams_mc_anticp, ith_cutoff)
    peptide_cv <- do_cv_peptides_mc(mer_cv)
    mer_res <- calc_cv_performance(mer_cv) %>% 
      mutate(layer = "mer")
    peptide_res <- calc_cv_performance(peptide_cv) %>% 
      mutate(layer = "peptide")
    print(paste0("done ", ith_cutoff))
    bind_rows(mer_res, peptide_res) %>% 
      mutate(cutoff = ith_cutoff)
  }) %>% bind_rows()
}
  

test_model_size <- function(cutoffs, mers, ngrams) {
  lapply(cutoffs, function(ith_cutoff) {
    imp_ngrams <- unique(unlist(unname(get_imp_ngrams_mc(ngrams_mc_anticp, mers_mc_anticp, ith_cutoff))))
    mer_model <- train_mc_model_mers(mers_mc_anticp, ngrams_mc_anticp, imp_ngrams)
    data.frame(n_imp_ngrams = length(imp_ngrams),
               model_size = print(format(object.size(mer_model), units = "MB", digits = 3)),
               stringsAsFactors = FALSE)
  }) %>% bind_rows()
}


test_cutoffs <- drake_plan(cv_res = test_cutoffs_cv(cutoffs),
                           model_sizes = test_model_size(cutoffs, mers_mc_anticp, ngrams_mc_anticp))


make(test_cutoffs, seed = 2938, cache = new_cache("testing_cutoffs"))
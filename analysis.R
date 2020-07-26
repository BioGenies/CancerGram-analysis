library(dplyr)
library(magrittr)
library(drake)
library(biogram)
library(seqinr)
library(pbapply)
library(ranger)
library(cvTools)
library(visNetwork)
library(hmeasure)
library(seqR)

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/Projekty/BioNgramProjects/CancerGram/"
}

source("./functions/raw_data.R")
source("./functions/cdhit_data.R")
source("./functions/nonstandard_AMPs.R")
source("./functions/cutting_seqs.R")
source("./functions/holdouts.R")
source("./functions/writing_benchmarks.R")
source("./functions/get_mers.R")
source("./functions/count_ampgrams.R")
source("./functions/do_cv.R")
source("./functions/train_model_peptides.R")

analysis_CancerGram <- drake_plan(gathered_data = gather_raw_data(),
                                  raw_data = read_raw_data(),
                                  cdhit_data = filter_cdhit(raw_data),
                                  negative_data = read_and_cut(data_path, lengths(cdhit_data)),
                                  cdhit_data_ids = generate_holdout_groups(cdhit_data),
                                  negative_data_ids = generate_holdout_groups(negative_data),
                                  benchmark_file = write_benchmark(pos = cdhit_data,
                                                                   pos_id = cdhit_data_ids,
                                                                   neg = negative_data,
                                                                   neg_id = negative_data_ids),
                                  mer_df = get_mers(pos = cdhit_data,
                                                    pos_id = cdhit_data_ids,
                                                    neg = negative_data,
                                                    neg_id = negative_data_ids),
                                  binary_ngrams = count_and_gather_ngrams(mer_df,
                                                                          c(1, rep(2, 4), rep(3, 4)),
                                                                          list(NULL, NULL, 1, 2, 3, c(0,0), c(0,1), c(1,0), c(1,1))),
                                  cv_raw = do_cv(mer_df, binary_ngrams))


make(analysis_CancerGram, seed = 2938)

library(dplyr)
library(magrittr)
library(drake)
library(biogram)
library(pbapply)
library(ranger)
library(cvTools)
library(visNetwork)
library(hmeasure)

if(Sys.info()[["nodename"]] %in% c("amyloid", "phobos", "huawei")) {
  data_path <- "/home/michal/Dropbox/AMP-analysis/AmpGram-analysis/"
}
if(Sys.info()[["nodename"]] %in% c("kasia-MACH-WX9", "ryzen")) {
  data_path <- "/home/kasia/Dropbox/AmpGram-analysis/"
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

analysis_CancerGram <- drake_plan(raw_data = read_raw_data(),
                                  cdhit_data = filter_cdhit(raw_data),
                                  negative_data = read_and_cut(data_path, lengths(cdhit_data)))


make(analysis_CancerGram, seed = 2938)

#' Processing sequence files
#' 
#' This function processess sequence files of AntiCP 2.0 datasets
#' and discards sequences shorter than 5 amino acids.
#' @param seq_file name of the file containing sequences
#' @return list of named sequences
process_sequences <- function(seq_file) {
  seqs <- readLines(paste0("./data/", seq_file))
  named_seqs <- seqs %>% 
    strsplit("") %>% 
    setNames(paste0(strsplit(seq_file, ".", fixed = TRUE)[[1]][1], "_", 1:length(seqs)))
  named_seqs[lengths(named_seqs) >= 5]
}

#' Get target
#' 
#' Adds target column to a data frame based on the peptide name.
#' @param df data frame to which target should be added
#' @return input data frame with additional column indicating the target name
get_target <- function(df) {
  mutate(df,
         target = factor(case_when(grepl("pos_train_main|pos_test_main|CancerPPD|AP|DRAMP", source_peptide) ~ "acp",
                                   grepl("neg_train_main|neg_test_main|dbAMP", source_peptide) ~ "amp",
                                   grepl("neg_train_alternate|neg_test_alternate", source_peptide) ~ "neg")))
}


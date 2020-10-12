# Processing AntiCP sequence files
process_sequences <- function(seq_file) {
  seqs <- readLines(paste0("./data/", seq_file))
  named_seqs <- seqs %>% 
    strsplit("") %>% 
    setNames(paste0(strsplit(seq_file, ".", fixed = TRUE)[[1]][1], "_", 1:length(seqs)))
  named_seqs[lengths(named_seqs) >= 5]
}

get_target <- function(df) {
  mutate(df,
         target = factor(case_when(grepl("pos_train_main|pos_test_main", source_peptide) ~ "acp",
                                   grepl("neg_train_main|neg_test_main", source_peptide) ~ "amp",
                                   grepl("neg_train_alternate|neg_test_alternate", source_peptide) ~ "neg")))
}


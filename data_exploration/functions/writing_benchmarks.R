#' Writes benchmark dataset into the file
#' 
#' Filters sequences from the positive and negative datasets for benchmark 
#' and writes them into the file in fasta format.
#' @param pos positive dataset 
#' @param pos_id IDs of the positive dataset
#' @param neg negative dataset
#' @param neg_id IDs of the negative dataset
#' @return a file containing sequences from benchmark dataset in fasta format
write_benchmark <- function(pos, pos_id, neg, neg_id) {
  seq_list <- c(pos[unlist(lapply(pos_id, function(ith_len_group) ith_len_group[["benchmark"]]))],
    neg[unlist(lapply(neg_id, function(ith_len_group) ith_len_group[["benchmark"]]))]) 
  write_fasta(seq_list, file = "results/benchmark.fasta") 
  print(paste0("Number of sequences in the benchmark dataset: ", length(seq_list))) 
}

# write_benchmark(pos = readd(cdhit_data),
#                 pos_id = readd(cdhit_data_ids),
#                 neg = readd(negative_data),
#                 neg_id = readd(negative_data_ids))

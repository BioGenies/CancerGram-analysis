library(dplyr)
library(ggplot2)
library(drake)

loadd(cdhit_data)
loadd(negative_data)

acp_lens <- data.frame(len = lengths(cdhit_data),
                       dataset = "Positive dataset")

ggplot(acp_lens, aes(x = dataset, y = len)) +
  geom_violin()

mutate(acp_lens, cutted_length = cut(len, breaks = c(0L:10 * 10, 150, 210))) %>% 
  ggplot(aes(x = cutted_length)) +
  geom_bar()


datasets <- list("positive" = cdhit_data, "negative" = negative_data)
aac <- lapply(1:length(datasets), function(i) {
  data.frame(table(unlist(datasets[[i]]))) %>% 
    setNames(c("aa", "freq")) %>% 
    mutate(freq = freq/sum(freq),
           dataset = factor(names(datasets[i]), levels = c("positive", "negative")))
}) %>% bind_rows()

ggplot(aac, aes(x = aa, y = freq, fill = dataset)) +
  geom_col(position = "dodge")

apd_raw <- read.csv("./data/apd_df.csv", stringsAsFactors = FALSE)

library(biogram)

AMP <- strsplit(apd_raw[!grepl("cancer", apd_raw[["Activity"]], ignore.case = TRUE), c("Sequence")], "")
nAMP_nANC <- biogram::read_fasta("./data/nonAMP_nonACP.fasta")
ANC <- cdhit_data

seq_list <- list(AMP = AMP, nAMP_nANC = nAMP_nANC, ANC = ANC)

prop_vec <- c(charge = "KLEP840101", 
              hydrophobicity = "KYTJ820101",
              hydrophobcity_black = "BLAS910101",
              positive = "FAUJ880111",
              helix_eq = "FINA770101",
              helix_prop = "KOEP990101")

encode_seq <- function(x, property) {
  sapply(x, function(ith_seq) {
    mean(aaprop[property, tolower(ith_seq)])
  })
}

prop_dat <- lapply(names(prop_vec), function(ith_prop_type) 
  lapply(names(seq_list), function(ith_seq_type) {
    data.frame(seq_type = ith_seq_type, prop_type = ith_prop_type, 
               value = encode_seq(seq_list[[ith_seq_type]], prop_vec[ith_prop_type]),
               stringsAsFactors = FALSE)
  }) %>% bind_rows
) %>% bind_rows()

ggplot(prop_dat, aes(x = seq_type, y = value)) +
  geom_boxplot() +
  facet_wrap(~prop_type, scales = "free_y")













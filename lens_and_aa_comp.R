library(dplyr)
library(ggplot2)
library(drake)
library(biogram)
library(tidyr)
library(seqinr)

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

#apd_raw <- read.csv("./data/apd_df.csv", stringsAsFactors = FALSE)

#AMP <- strsplit(apd_raw[!grepl("cancer", apd_raw[["Activity"]], ignore.case = TRUE), c("Sequence")], "")

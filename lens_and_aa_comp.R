library(dplyr)
library(ggplot2)
library(drake)

loadd(cdhit_data)
loadd(negative_data)

acp_lens <- data.frame(len = lengths(cdhit_data),
                       dataset = "Positive dataset")

ggplot(acp_lens, aes(x = dataset, y = len)) +
  geom_violin()

datasets <- list("positive" = cdhit_data, "negative" = negative_data)
aac <- lapply(1:length(datasets), function(i) {
  data.frame(table(unlist(datasets[[i]]))) %>% 
    setNames(c("aa", "freq")) %>% 
    mutate(freq = freq/sum(freq),
           dataset = factor(names(datasets[i]), levels = c("positive", "negative")))
}) %>% bind_rows()


ggplot(aac, aes(x = aa, y = freq, fill = dataset)) +
  geom_col(position = "dodge")
  
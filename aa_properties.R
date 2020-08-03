library(dplyr)
library(ggplot2)
library(drake)
library(biogram)
library(tidyr)
library(seqinr)


ANC <- readd(cdhit_data)
AmpGram_data <- readd(cdhit_data, path = "/home/kasia/RProjects/AmpGram-analysis/.drake")
AMP <- AmpGram_data[which(!(AmpGram_data %in% ANC))]
nAMP_nANC <- readd(negative_data)

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

filter(prop_dat, prop_type %in% c("charge", "hydrophobicity")) %>% 
  group_by(prop_type) %>% 
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = prop_type, values_from = value) %>% 
  ggplot(aes(x = charge, y = hydrophobicity, color = seq_type, fill = seq_type)) +
  #  geom_density2d(aes(alpha = ..level..))
  stat_density_2d(aes(alpha = ..level..), geom = "polygon", color = "black", size = 0.4) +
  facet_wrap(~seq_type)



### Testing all properties on both datasets
all_seq_list <- list(AMP = readd(amp), neg = readd(negative_data), ACP = readd(cdhit_data), 
                     AMP_AntiCP = readd(neg_train_main), neg_AntiCP = readd(neg_train_alt), 
                     ACP_AntiCP = readd(pos_train_main))
data(aaindex)

all_prop_dat <- lapply(names(aaindex), function(ith_prop_type) 
  lapply(names(all_seq_list), function(ith_seq_type) {
    data.frame(seq_type = gsub("_AntiCP", "", ith_seq_type), prop_type = ith_prop_type, 
               value = encode_seq(all_seq_list[[ith_seq_type]], ith_prop_type),
               dataset = ifelse(grepl("AntiCP", ith_seq_type), "AntiCP", "CancerGram"),
               stringsAsFactors = FALSE)
  }) %>% bind_rows
) %>% bind_rows()

diffs <- group_by(all_prop_dat, dataset, prop_type, seq_type) %>% 
  summarize(mean_val = mean(value)) %>% 
  pivot_wider(names_from = seq_type, values_from = mean_val) %>% 
  group_by(prop_type, dataset) %>%
  summarise(diff_acp_amp = abs(ACP - AMP)) %>% 
  arrange(desc(diff_acp_amp)) 

prop_order <- unique(diffs[["prop_type"]])
dataset_color = c(ACP = "#d73027", AMP = "#fc8d59", neg = "#91bfdb")

lapply(seq(1, length(prop_order), 8), function(i) {
  props <- prop_order[i:(i+7)]
  plot_dat <- filter(all_prop_dat, prop_type %in% props) %>% 
    group_by(dataset) 
  plot <- ggplot(plot_dat, aes(x = seq_type, y = value, fill = seq_type)) +
    geom_violin() +
    scale_fill_manual("Dataset", values = dataset_color) +
    facet_grid(dataset ~ prop_type, 
               labeller = labeller(prop_type = sapply(props, function(ith_prop) gsub("(", "\n(", aaindex[[ith_prop]][["D"]], fixed = TRUE))))
  ggsave(paste0("./prop_plots/props_", i, "-", i+7, ".png"), plot, "png", width = 28, height = 10)
})



### Check properties of incorrectly predicted peptides

# Properties for whole datasets
calc_prop_values <- function(seq_list, prop_vec) {
  lapply(names(prop_vec), function(ith_prop_type) {
    lapply(names(seq_list), function(ith_seq_type) {
      data.frame(seq_type = ith_seq_type, prop_type = ith_prop_type, 
                 value = encode_seq(seq_list[[ith_seq_type]], prop_vec[ith_prop_type]),
                 stringsAsFactors = FALSE)
    }) %>% bind_rows
  }) %>% bind_rows()
}

# Properties for selected peptides
get_props_for_selected <- function(selected_preds, prop_vec) {
  lapply(selected_preds, function(ith_preds) {
    seqs <- c(pos_test_main, neg_test_main, neg_test_alt)[ith_preds[["source_peptide"]]]
    lapply(names(prop_vec), function(ith_prop_type) {
      lapply(names(seqs), function(ith_seq) {
        data.frame(seq = ith_seq, prop_type = ith_prop_type, 
                   value = encode_seq(seqs[ith_seq], prop_vec[ith_prop_type]),
                   stringsAsFactors = FALSE)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
}

pred_res <- readd(benchmark_peptide_preds_mc_anticp)

acp_as_neg <- filter(pred_res, target == "acp" & neg > 0.5) #14
neg_as_acp <- filter(pred_res, target == "neg" & acp > 0.5) #1
amp_as_acp <- filter(pred_res, target == "amp" & acp > 0.5) #20
acp_as_amp <- filter(pred_res, target == "acp" & amp > 0.5) #40

seq_list <- list(AMP = readd(pos_train_main), neg = readd(neg_train_alt), ACP = readd(neg_train_main)) 

prop_dat <- calc_prop_values(seq_list, prop_vec)

acp_neg <- get_props_for_selected(list(acp_as_neg = acp_as_neg, neg_as_acp = neg_as_acp), prop_vec) %>% 
  mutate(seq_type = case_when(grepl("pos_test_main", seq) ~ "ACP",
                              grepl("neg_test_main", seq) ~ "AMP",
                              grepl("neg_test_alt", seq) ~ "neg"))

acp_amp <- get_props_for_selected(list(acp_as_amp = acp_as_amp, amp_as_acp = amp_as_acp), prop_vec) %>% 
  mutate(seq_type = case_when(grepl("pos_test_main", seq) ~ "ACP",
                              grepl("neg_test_main", seq) ~ "AMP",
                              grepl("neg_test_alt", seq) ~ "neg"))

ggplot(prop_dat, aes(x = seq_type, y = value, fill = seq_type)) +
  geom_violin() +
  scale_fill_manual("Dataset", values = dataset_color) +
  facet_wrap(~prop_type, scales = "free_y") +
  #geom_point(data = acp_neg, position = "jitter")
  geom_point(data = acp_amp, position = "jitter")


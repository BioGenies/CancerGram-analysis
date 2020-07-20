library(dplyr)
library(ggplot2)
library(drake)

source("functions/do_cv.R")

cv_raw <- readRDS("cv_raw.RDS")

perf_res <- lapply(unique(cv_raw[["comb"]]), function(ith_comb) 
  lapply(1L:5, function(ith_fold) 
    lapply(unique(cv_raw[["group"]]), function(ith_group) {
      perf_dat <- filter(cv_raw, 
                         group == ith_group, 
                         fold == ith_fold,
                         comb == ith_comb) %>% 
        mutate(Decision = factor(ifelse(pred >= 0.5, "TRUE", "FALSE")),
               target = as.factor(target))

      data.frame(AUC = mlr3measures::auc(perf_dat[["target"]], perf_dat[["pred"]], "TRUE"),
                 MCC = mlr3measures::mcc(perf_dat[["target"]], perf_dat[["Decision"]], "TRUE"),
                 Precision = mlr3measures::precision(perf_dat[["target"]], perf_dat[["Decision"]], "TRUE"),
                 Sensitivity = mlr3measures::sensitivity(perf_dat[["target"]], perf_dat[["Decision"]], "TRUE"),
                 Specificity = mlr3measures::specificity(perf_dat[["target"]], perf_dat[["Decision"]], "TRUE")) %>% 
        mutate(fold = ith_fold, comb = length(strsplit(ith_comb, ",")[[1]])/2, group = ith_group)
    }) %>% bind_rows
  ) %>% bind_rows
) %>% 
  bind_rows %>% 
  mutate(group = factor(group, levels = sort_group(unique(group))))


ggplot(perf_res, aes(x = group, y = AUC)) +
  geom_point() +
  stat_summary(fun = mean, geom = "point", size = 4, color = "red") +
  facet_wrap(~ comb, nrow = 1)

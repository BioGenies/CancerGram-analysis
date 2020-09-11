get_cv_pred_table <- function(cv_performance, outfile_name) {
  tab <- cv_performance %>% 
    pivot_longer(AUC:Accuracy, names_to = "Measure", values_to = "Value") %>% 
    group_by(dataset, Measure) %>%
    summarise(mean = mean(Value),
              sd = sd(Value)) %>% 
    mutate(value = paste0(round(mean, 2), " (+/-", round(sd, 3), ")")) %>% 
    select(-c(mean, sd)) %>%
    pivot_wider(names_from = dataset, values_from = value)
  write.csv(tab, paste0(outfile_name, ".csv"), row.names = FALSE)
    xtable(tab, caption = "", label = "", align = "cccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
      writeLines(paste0(outfile_name, ".txt"))
}

get_datasets_table <- function(acp_train, acp_test, amp_train, amp_test, neg_train, neg_test) {
  df <- data.frame(Dataset = c("Train-test", "Benchmark"),
             ACP = c(length(acp_train), length(acp_test)),
             AMP = c(length(amp_train), length(amp_test)),
             neg = c(length(neg_train), length(neg_test)))
  write.csv(df, paste0(data_path, "datasets_table.csv"), row.names = FALSE)
  xtable(df, caption = "", label = "", align = "ccccc") %>% 
    print(include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top", label.placement = "top") %>% 
    writeLines(paste0(data_path, "datasets_table.txt"))
}
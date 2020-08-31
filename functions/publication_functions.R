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
    print(include.rownames = FALSE) %>% 
      writeLines(paste0(outfile_name, ".txt"))
}


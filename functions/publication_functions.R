get_cv_pred_table <- function(cv_performance) {
  cv_performance %>% 
    pivot_longer(AUC:Accuracy, names_to = "Measure", values_to = "Value") %>% 
    group_by(dataset, Measure) %>%
    summarise(mean = mean(Value),
              sd = sd(Value)) %>% 
    mutate(value = paste0(round(mean, 3), " (+/-", round(sd, 3), ")")) %>% 
    select(-c(mean, sd)) %>% 
    pivot_wider(names_from = Measure, values_from = value) %>% 
    xtable()

}
library(dplyr)

model1 <- read.csv("./data/anticp_preds_model1.csv", stringsAsFactors = FALSE) %>% 
  mutate(target = "neg",
         model = "ACP/AMP")
model2 <- read.csv("./data/anticp_preds_model2.csv", stringsAsFactors = FALSE) %>% 
  mutate(target = "amp",
         model = "ACP/neg")

bind_rows(model1, model2) %>% 
  group_by(target, model, Prediction) %>% 
  summarise(count = n())
  
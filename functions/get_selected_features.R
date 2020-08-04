get_selected_features <- function(alphabets) {
  selected_features <- alphabets[["features"]]
  data("aaindex")
  prop_names <- sapply(1:length(aaindex), function(x) {
    paste(attributes(aaindex[x]))
  })
  properties <- lapply(prop_names, function(x){
    paste(aaindex[[x]][["D"]])
  }) %>% unlist() 
  aaprop_table <- data.frame(name = prop_names, feature = properties, stringsAsFactors = FALSE) %>% 
    mutate(selected = ifelse(name %in% selected_features, "yes", ""))
  write.csv(aaprop_table, "reports/selected_features.csv", row.names = FALSE)
}

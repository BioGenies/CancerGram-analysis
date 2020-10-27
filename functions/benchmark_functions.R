get_independent_dataset_preds <- function(dataset, anticp_preds_file) {

anticp_res <- read.csv(anticp_preds_file, stringsAsFactors = FALSE)  %>% 
  mutate(Decision = factor(ifelse(Prediction == "AntiCP", "acp", "amp")),
         Software = "AntiCP 2.0") %>%
  select(c(`X..Sequence_ID`, Score, Decision, Software)) %>% 
  setNames(c("source_peptide", "Prob", "Decision", "Software")) %>% 
  get_target() %>% 
  mutate(target = factor(ifelse(target == "acp", "acp", "amp")))

validation_mers <- mer_df_from_list(dataset)
validation_ngrams <- count_imp_ngrams(validation_mers, imp_ngrams)
validation_mer_preds <- cbind(validation_mers, predict(mer_model, 
                                                       as.matrix(validation_ngrams))[["predictions"]])
validation_stats <- mutate(calculate_statistics_mc(validation_mer_preds, c("acp", "amp", "neg")),
                           target = factor(case_when(grepl("CUTTED", source_peptide) ~ "neg",
                                                     grepl("dbAMP", source_peptide) ~ "amp",
                                                     grepl("CancerPPD|AP|DRAMP", source_peptide) ~ "acp")))
validation_peptide_preds <- cbind(validation_stats[, c("source_peptide", "target")],
                                  predict(peptide_model, validation_stats)[["predictions"]])

cancergram_res <- validation_peptide_preds %>% 
  mutate(Decision = factor(case_when(acp > amp & acp > neg  ~ "acp",
                                     amp > acp & amp > neg ~ "amp",
                                     neg > amp & neg > acp ~ "amp")),
         Software = "CancerGram",
         target = factor(ifelse(target == "acp", "acp", "amp"))) %>% 
  select(c(source_peptide, acp, target, Decision, Software)) %>% 
  setNames(c("source_peptide", "Prob", "target", "Decision", "Software")) 

get_metrics(bind_rows(anticp_res, cancergram_res))

}


get_metrics <- function(preds) {
  preds %>% 
    group_by(Software) %>% 
    summarise(MCC = mlr3measures::mcc(target, Decision, "acp"),
              Precision = mlr3measures::precision(target, Decision, "acp"),
              Sensitivity = mlr3measures::sensitivity(target, Decision, "acp"),
              Specificity = mlr3measures::specificity(target, Decision, "acp"),
              Accuracy = mlr3measures::acc(target, Decision),
              AUC = mlr3measures::auc(target, Prob, "acp"))
}



get_validation_mACPpred_dataset_preds <- function(pos_mACPpred_dataset, neg_mACPpred_dataset,
                                                  pos_validation, neg_validation, 
                                                  cancergram_preds, mACPpred_preds) {
  macp_p <- read_fasta(pos_mACPpred_dataset)
  macp_n <- read_fasta(neg_mACPpred_dataset)
  pos_test_main[which(!(pos_test_main %in% macp_p))] %>% 
    write_fasta("./data/pos_test_main_without_mACPpred.fa")
  macp_duplicates <- names(pos_test_main[which(pos_test_main %in% macp_p)])
  cancergram_preds <- filter(benchmark_peptide_preds, !(source_peptide %in% macp_duplicates), target %in% c("acp", "amp")) %>% 
    mutate(Decision = factor(case_when(acp > amp & acp > neg  ~ "acp",
                                       amp > acp & amp > neg ~ "amp",
                                       neg > amp & neg > acp ~ "amp")),
           Software = "CancerGram",
           target = factor(ifelse(target == "acp", "acp", "amp"))) %>% 
    select(source_peptide, acp, target, Decision, Software) %>% 
    setNames(c("source_peptide", "Prob", "target", "Decision", "Software"))
  
  mACP_preds <- read.csv(mACPpred_preds, stringsAsFactors = FALSE) %>% 
    mutate(source_peptide = gsub("-", "_", `FASTA.ID`),
           Decision = factor(ifelse(`ACP.or.Non.ACP` == "ACP", "acp", "amp")),
           Software = "mACPpred") %>% 
    select(source_peptide, `Prob`, Decision, Software) %>% 
    get_target() %>% 
    mutate(target = factor(ifelse(target == "acp", "acp", "amp")))
  
  get_metrics(bind_rows(cancergram_preds, mACP_preds))
  
}
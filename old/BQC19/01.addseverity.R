setwd("/project/richards/restricted/bqc19_release6/Clinical_Phenotypic")

classfication <- fread("bqc_blood_samples_with_Inf_nonInf_classification.tsv")
pheno_raw <- read.csv("redcap_clinical_data_raw_2022-03-24.csv", header = T, fileEncoding="latin1", sep="\t")

###severity final
icu_id <- pheno_raw %>% filter(icu == 1) %>% dplyr::select(BQCID) %>% unique()
hosp_id <- pheno_raw %>% filter(admit == 1) %>% dplyr::select(BQCID) %>% unique()
death_id <- pheno_raw %>% filter(death_dc == 1 | death == 1 | exit_reason == 0 | visit_status_spec == 3 |
                                 interview_whynot == 0) %>% dplyr::select(BQCID) %>% unique()
resp_id <- pheno_raw %>% filter(respsupport___2 == 1 | respsupport___3 == 1 | respsupport___4 == 1 |
                                  othersupport___4 == 1 | c_ards == 1 | covid_severity %in% c(3,4)) %>%
  dplyr::select(BQCID) %>% unique()
  
classfication <- classfication %>% mutate(A2 = case_when(covid19_test_result == "Positive" & BQCID %in% icu_id$BQCID ~ 1,
                                                         covid19_test_result == "Positive" & BQCID %in% death_id$BQCID ~ 1,
                                                         covid19_test_result == "Positive" & BQCID %in% resp_id$BQCID ~ 1,
                                                         TRUE ~ 0),
                                          B2 = case_when(A2 == 1 ~ 1,
                                                         covid19_test_result == "Positive" & BQCID %in% hosp_id$BQCID ~ 1,
                                                         TRUE ~ 0),
                                          C2 = case_when(A2 == 1 ~ 1,
                                                         B2 == 1 ~ 1,
                                                         covid19_test_result == "Positive" ~ 1,
                                                         covidstatus == "Positive" ~ 1,
                                                         TRUE ~ 0))
classfication <- classfication %>% mutate(final_severity = case_when(covid19_test_result == "Positive" & BQCID %in% icu_id$BQCID ~ "SEVERE",
                                                                     covid19_test_result == "Positive" & BQCID %in% death_id$BQCID ~ "SEVERE",
                                                                     covid19_test_result == "Positive" & BQCID %in% resp_id$BQCID ~ "SEVERE",
                                                                     covid19_test_result == "Positive" & BQCID %in% hosp_id$BQCID ~ "MODERATE",
                                                                     covid19_test_result == "Positive" ~ "MILD",
                                                                     covidstatus == "Negative" ~ "0NonCOVID"))

write.table(classfication, file="/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/05.BQC/phenotype/release6/bqc_blood_samples_with_Inf_nonInf_classification_with_severity.tsv", quote=F, col.names = T, row.names = F, sep="\t")

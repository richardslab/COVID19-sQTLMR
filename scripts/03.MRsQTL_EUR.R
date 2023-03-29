setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

exp_sQTL_Lung_EUR <- readRDS("exposure_sQTL_Lung_EUR.rds")
exp_sQTL_WBC_EUR <- readRDS("exposure_sQTL_WBC_EUR.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/A2_EUR.r7.b38.rds")

dat <- harmonise_data(
  exposure_dat = exp_sQTL_Lung_EUR, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))
dat %>% saveRDS("dat_sQTL_Lung_A2_EUR.rds")
res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

saveRDS(res, "A2_Lung_sQTL_EUR.rds")

dat <- harmonise_data(
  exposure_dat = exp_sQTL_WBC_EUR, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))
dat %>% saveRDS("dat_sQTL_WBC_A2_EUR.rds")

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

saveRDS(res, "A2_WBC_sQTL_EUR.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/B2_EUR.r7.b38.rds")
dat <- harmonise_data(
  exposure_dat = exp_sQTL_Lung_EUR, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))
dat %>% saveRDS("dat_sQTL_Lung_B2_EUR.rds")

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

saveRDS(res, "B2_Lung_sQTL_EUR.rds")

dat <- harmonise_data(
  exposure_dat = exp_sQTL_WBC_EUR, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))
dat %>% saveRDS("dat_sQTL_WBC_B2_EUR.rds")


res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

saveRDS(res, "B2_WBC_sQTL_EUR.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/C2_EUR.r7.b38.rds")

#OUTCOME
dat <- harmonise_data(
  exposure_dat = exp_sQTL_Lung_EUR, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))
dat %>% saveRDS("dat_sQTL_Lung_C2_EUR.rds")

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))


saveRDS(res, "C2_Lung_sQTL_EUR.rds")

dat <- harmonise_data(
  exposure_dat = exp_sQTL_WBC_EUR, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))
dat %>% saveRDS("dat_sQTL_WBC_C2_EUR.rds")


res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

saveRDS(res, "C2_WBC_sQTL_EUR.rds")

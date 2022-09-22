setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

#sQTL
sQTL_WBC <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Whole_Blood.v8.independent_sqtls.txt.gz")
sQTL_WBC <- sQTL_WBC %>% arrange(pval_nominal)
sQTL_WBC <- sQTL_WBC %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
sQTL_WBC <- sQTL_WBC %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                N = round(ma_count/maf*(1/2)))

sQTL_WBC <- sQTL_WBC %>% mutate(CHRPOS = gsub("chrX", "chr23", CHRPOS))

exp_sQTL_WBC <- format_data(sQTL_WBC, type="exposure",
                       phenotype_col = "phenotype_id",
                       snp_col = "CHRPOS",
                       beta_col = "slope",
                       se_col = "slope_se",
                       effect_allele_col = "EA",
                       other_allele_col = "NEA",
                       pval_col = "pval_nominal",
                       samplesize_col = "N",
                       eaf = "maf",
                       chr_col = "CHR",
                       pos_col = "POS"
)

out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/C2_EUR.r7.b38.rds")

tmp <- out %>% select(SNP, eaf.outcome, effect_allele.outcome, other_allele.outcome)
exp_sQTL_WBC1 <- exp_sQTL_WBC %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_sQTL_WBC1 <- exp_sQTL_WBC1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                   effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                   TRUE ~ eaf.exposure))
exp_sQTL_WBC1 <- exp_sQTL_WBC1 %>% select(colnames(exp_sQTL_WBC))


##lung
sQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Lung.v8.independent_sqtls.txt.gz")
sQTL_Lung <- sQTL_Lung %>% arrange(pval_nominal)
sQTL_Lung <- sQTL_Lung %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                  POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                  EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                  NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
sQTL_Lung <- sQTL_Lung %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                  N = round(ma_count/maf*(1/2)))

exp_sQTL_Lung <- format_data(sQTL_Lung, type="exposure",
                       phenotype_col = "phenotype_id",
                       snp_col = "CHRPOS",
                       beta_col = "slope",
                       se_col = "slope_se",
                       effect_allele_col = "EA",
                       other_allele_col = "NEA",
                       pval_col = "pval_nominal",
                       samplesize_col = "N",
                       eaf = "maf",
                       chr_col = "CHR",
                       pos_col = "POS"
)

exp_sQTL_Lung1 <- exp_sQTL_Lung %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_sQTL_Lung1 <- exp_sQTL_Lung1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                   effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                   TRUE ~ eaf.exposure))
exp_sQTL_Lung1 <- exp_sQTL_Lung1 %>% select(colnames(exp_sQTL_Lung))


#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/A2_EUR.r7.b38.rds")

dat <- harmonise_data(
  exposure_dat = exp_sQTL_WBC1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "A2_WBC_sQTL.rds")

dat <- harmonise_data(
  exposure_dat = exp_sQTL_Lung1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "A2_Lung_sQTL.rds")


#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/B2_EUR.r7.b38.rds")
dat <- harmonise_data(
  exposure_dat = exp_sQTL_WBC1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "B2_WBC_sQTL.rds")


#OUTCOME
dat <- harmonise_data(
  exposure_dat = exp_sQTL_Lung1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "B2_Lung_sQTL.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/C2_EUR.r7.b38.rds")
dat <- harmonise_data(
  exposure_dat = exp_sQTL_WBC1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "C2_WBC_sQTL.rds")


#OUTCOME
dat <- harmonise_data(
  exposure_dat = exp_sQTL_Lung1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "C2_Lung_sQTL.rds")


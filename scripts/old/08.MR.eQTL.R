setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

#eQTL
eQTL_WBC <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_Analysis_v8_eQTL_independent/Whole_Blood.v8.independent_eqtls.txt.gz")
eQTL_WBC <- eQTL_WBC %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
eQTL_WBC <- eQTL_WBC %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                N = round(minor_allele_count/maf*(1/2)))

eQTL_WBC <- eQTL_WBC %>% mutate(CHRPOS = gsub("chrX", "chr23", CHRPOS))
eQTL_WBC <- eQTL_WBC %>% filter(pval_nominal < 1e-5)

exp_eQTL_WBC <- format_data(eQTL_WBC, type="exposure",
                            phenotype_col = "gene_id",
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
saveRDS(exp_eQTL_WBC, file="exposure_eQTL_WBC.rds")


eQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_Analysis_v8_eQTL_independent/Lung.v8.independent_eqtls.txt.gz")
eQTL_Lung <- eQTL_Lung %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
eQTL_Lung <- eQTL_Lung %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                N = round(minor_allele_count/maf*(1/2)))

eQTL_Lung <- eQTL_Lung %>% mutate(CHRPOS = gsub("chrX", "chr23", CHRPOS))
eQTL_Lung <- eQTL_Lung %>% filter(pval_nominal < 1e-5)


exp_eQTL_Lung <- format_data(eQTL_Lung, type="exposure",
                            phenotype_col = "gene_id",
                            snp_col = "CHRPOS",
                            beta_col = "slope",
                            se_col = "slope_se",
                            effect_allele_col = "EA",
                            other_allele_col = "NEA",
                            eaf = "maf",
                            pval_col = "pval_nominal",
                            samplesize_col = "N",
                            chr_col = "CHR",
                            pos_col = "POS",
)

saveRDS(exp_eQTL_Lung, file="exposure_eQTL_Lung.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/A2_ALL.r7.b38.rds")
tmp <- out %>% select(SNP, eaf.outcome, effect_allele.outcome, other_allele.outcome)
exp_eQTL_WBC1 <- exp_eQTL_WBC %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_eQTL_WBC1 <- exp_eQTL_WBC1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                   effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                   TRUE ~ eaf.exposure))
exp_eQTL_WBC1 <- exp_eQTL_WBC1 %>% select(colnames(exp_eQTL_WBC))
dat <- harmonise_data(
  exposure_dat = exp_eQTL_WBC1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res %>% arrange(p) %>% head()

saveRDS(res, "A2_WBC_eQTL.rds")

exp_eQTL_Lung1 <- exp_eQTL_Lung %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_eQTL_Lung1 <- exp_eQTL_Lung1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                   effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                   TRUE ~ eaf.exposure))
exp_eQTL_Lung1 <- exp_eQTL_Lung1 %>% select(colnames(exp_eQTL_Lung))
dat <- harmonise_data(
  exposure_dat = exp_eQTL_Lung1, 
  outcome_dat = out, action = 2
)


res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "A2_Lung_eQTL.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/B2_ALL.r7.b38.rds")
tmp <- out %>% select(SNP, eaf.outcome, effect_allele.outcome, other_allele.outcome)
exp_eQTL_WBC1 <- exp_eQTL_WBC %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_eQTL_WBC1 <- exp_eQTL_WBC1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                   effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                   TRUE ~ eaf.exposure))
exp_eQTL_WBC1 <- exp_eQTL_WBC1 %>% select(colnames(exp_eQTL_WBC))
dat <- harmonise_data(
  exposure_dat = exp_eQTL_WBC1, 
  outcome_dat = out, action = 2
)

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "B2_WBC_eQTL.rds")


#OUTCOME
exp_eQTL_Lung1 <- exp_eQTL_Lung %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_eQTL_Lung1 <- exp_eQTL_Lung1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                     effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                     TRUE ~ eaf.exposure))
exp_eQTL_Lung1 <- exp_eQTL_Lung1 %>% select(colnames(exp_eQTL_Lung))
dat <- harmonise_data(
  exposure_dat = exp_eQTL_Lung1, 
  outcome_dat = out, action = 2
)


res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "B2_Lung_eQTL.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/C2_ALL.r7.b38.rds")
tmp <- out %>% select(SNP, eaf.outcome, effect_allele.outcome, other_allele.outcome)
exp_eQTL_WBC1 <- exp_eQTL_WBC %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_eQTL_WBC1 <- exp_eQTL_WBC1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                   effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                   TRUE ~ eaf.exposure))
exp_eQTL_WBC1 <- exp_eQTL_WBC1 %>% select(colnames(exp_eQTL_WBC))
dat <- harmonise_data(
  exposure_dat = exp_eQTL_WBC1, 
  outcome_dat = out, action = 2
)


res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "C2_WBC_eQTL.rds")


#OUTCOME
exp_eQTL_Lung1 <- exp_eQTL_Lung %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_eQTL_Lung1 <- exp_eQTL_Lung1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                     effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                     TRUE ~ eaf.exposure))
exp_eQTL_Lung1 <- exp_eQTL_Lung1 %>% select(colnames(exp_eQTL_Lung))
dat <- harmonise_data(
  exposure_dat = exp_eQTL_Lung1, 
  outcome_dat = out, action = 2
)
res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "C2_Lung_eQTL.rds")


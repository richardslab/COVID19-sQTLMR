setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)


#OUTCOME
a2_b38 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/main/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz")
a2_b38 <- a2_b38 %>% mutate(SNP_POS = paste0("chr",`#CHR`,":",POS))
out_dat <- format_data(a2_b38, type="outcome",
                       snp_col = "SNP_POS",
                       effect_allele_col = "ALT",
                       other_allele_col = "REF",
                       eaf_col = "all_meta_AF",
                       beta_col = "all_inv_var_meta_beta",
                       se_col = "all_inv_var_meta_sebeta",
                       pval_col = "all_inv_var_meta_p")

saveRDS(out_dat, file="/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/A2_ALL.r7.b38.rds")

b2_b38 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/main/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz")
b2_b38 <- b2_b38 %>% mutate(SNP_POS = paste0("chr",`#CHR`,":",POS))
out_dat <- format_data(b2_b38, type="outcome",
                       snp_col = "SNP_POS",
                       effect_allele_col = "ALT",
                       other_allele_col = "REF",
                       eaf_col = "all_meta_AF",
                       beta_col = "all_inv_var_meta_beta",
                       se_col = "all_inv_var_meta_sebeta",
                       pval_col = "all_inv_var_meta_p")

saveRDS(out_dat, file="/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/B2_ALL.r7.b38.rds")

b1_b38 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/main/COVID19_HGI_B1_ALL_leave_23andme_20220403.tsv.gz")
b1_b38 <- b1_b38 %>% mutate(SNP_POS = paste0("chr",`#CHR`,":",POS))
out_dat <- format_data(b1_b38, type="outcome",
                       snp_col = "SNP_POS",
                       effect_allele_col = "ALT",
                       other_allele_col = "REF",
                       eaf_col = "all_meta_AF",
                       beta_col = "all_inv_var_meta_beta",
                       se_col = "all_inv_var_meta_sebeta",
                       pval_col = "all_inv_var_meta_p")

saveRDS(out_dat, file="/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/B1_ALL.r7.b38.rds")


c2_b38 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/main/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz")
c2_b38 <- c2_b38 %>% mutate(SNP_POS = paste0("chr",`#CHR`,":",POS))
out_dat <- format_data(c2_b38, type="outcome",
                       snp_col = "SNP_POS",
                       effect_allele_col = "ALT",
                       other_allele_col = "REF",
                       eaf_col = "all_meta_AF",
                       beta_col = "all_inv_var_meta_beta",
                       se_col = "all_inv_var_meta_sebeta",
                       pval_col = "all_inv_var_meta_p")

saveRDS(out_dat, file="/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/C2_ALL.r7.b38.rds")

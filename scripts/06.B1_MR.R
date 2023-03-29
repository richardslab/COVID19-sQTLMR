setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

exp_sQTL_WBC <- readRDS("exposure_sQTL_WBC_EUR.rds")
exp_sQTL_Lung <- readRDS("exposure_sQTL_Lung_EUR.rds")

#OUTCOME
out <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/B1_ALL.r7.b38.rds")
tmp <- out %>% select(SNP, eaf.outcome, effect_allele.outcome, other_allele.outcome)
exp_sQTL_WBC1 <- exp_sQTL_WBC %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_sQTL_WBC1 <- exp_sQTL_WBC1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                   effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                   TRUE ~ eaf.exposure))
exp_sQTL_WBC1 <- exp_sQTL_WBC1 %>% select(colnames(exp_sQTL_WBC))
dat <- harmonise_data(
  exposure_dat = exp_sQTL_WBC1, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res %>% arrange(p) %>% head()

saveRDS(res, "B1_WBC_sQTL_ALL.rds")

exp_sQTL_Lung1 <- exp_sQTL_Lung %>% inner_join(tmp, by=c("SNP"="SNP"))
exp_sQTL_Lung1 <- exp_sQTL_Lung1 %>% mutate(eaf.exposure = case_when(effect_allele.exposure == effect_allele.outcome ~ eaf.exposure,
                                                                     effect_allele.exposure == other_allele.outcome ~ 1 - eaf.exposure,
                                                                     TRUE ~ eaf.exposure))
exp_sQTL_Lung1 <- exp_sQTL_Lung1 %>% select(colnames(exp_sQTL_Lung))
dat <- harmonise_data(
  exposure_dat = exp_sQTL_Lung1, 
  outcome_dat = out, action = 2
)

dat <- dat %>% mutate(SNP = paste0(SNP,":",effect_allele.exposure,":",other_allele.exposure))

res <- mr_singlesnp(dat,
                    single_method = "mr_wald_ratio",
                    all_method = c("mr_ivw", "mr_egger_regression"))

res <- res %>% mutate(p.adj = p*length(unique(dat$exposure)))
res %>% arrange(p) %>% filter(p.adj < 0.05)

saveRDS(res, "B1_Lung_sQTL_ALL.rds")


B1_Lung_sQTL <- readRDS("B1_Lung_sQTL_ALL.rds")
B1_Lung_sQTL <- B1_Lung_sQTL %>% drop_na(p) %>% mutate(Tissue = "Lung")

B1_WBC_sQTL <- readRDS("B1_WBC_sQTL_ALL.rds")
B1_WBC_sQTL <- B1_WBC_sQTL %>% drop_na(p) %>% mutate(Tissue = "Whole Blood")
tmp <- readRDS("Whole_Blood_cluster_intron_gene_map.rds")
tmp <- tmp %>% mutate(exposure = paste0(pos,":",gene),
                      exposure1 = paste0(pos,":",cluster,":",gene)) %>% dplyr::select(exposure, exposure1)
B1_WBC_sQTL <- B1_WBC_sQTL %>% inner_join(tmp, by="exposure")
B1_WBC_sQTL <- B1_WBC_sQTL %>% mutate(exposure = exposure1)
B1_WBC_sQTL %>% saveRDS("B1_WBC_sQTL_ALL.rds")

ALL <- bind_rows(B1_WBC_sQTL, B1_Lung_sQTL)
ALL <- ALL %>% mutate(OR = exp(b),
         LL = exp(b + qnorm(0.025)*se),
         UL = exp(b + qnorm(0.975)*se),
         ensembleID = str_split(exposure, pattern="\\:", simplify = T)[,5],
         ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])

ALL_Lung_sQTL_EUR <- readRDS("EUR_Lung_sQTL_EUR_coloc.rds")
ALL_WBC_sQTL_EUR <- readRDS("EUR_WBC_sQTL_EUR_coloc.rds")

tmp_lung <- ALL_Lung_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8)
tmp_lung %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_lung %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_lung %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()

tmp_wbc <- ALL_WBC_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8)
tmp_wbc %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_wbc %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_wbc %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()


tmp <- bind_rows(tmp_lung, tmp_wbc)


ALL <- ALL %>% filter(ensembleID %in% c("ENSG00000068650", "ENSG00000142002", "ENSG00000089127", "ENSG00000168743", "ENSG00000160783"))
library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(ALL$ensembleID))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

ALL <- ALL %>% left_join(G_list, by=c("ensembleID"="ensembl_gene_id"))

final <- ALL %>% dplyr::select(exposure, ensembleID, hgnc_symbol, Tissue, SNP, OR, LL, UL, p)

final <- final %>% filter(!grepl("chr19:4714337:4719851", exposure))

final <- final %>% filter(!grepl("chr13:112875941:112881816", exposure))


final <- final %>% filter(!grepl("chr4:105959126:105968895", exposure))

final <- final %>% filter(!grepl("chr1:156233728:156236349", exposure))


final %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable3.B1.xlsx")



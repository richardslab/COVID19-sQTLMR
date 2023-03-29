setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

#eQTL
res_a <- readRDS("A2_Lung_eQTL_EUR.rds") %>% mutate(outcome = "A2")
res_b2 <- readRDS("B2_Lung_eQTL_EUR.rds") %>% mutate(outcome = "B2")
res_c <- readRDS("C2_Lung_eQTL_EUR.rds") %>% mutate(outcome = "C2")

res <- bind_rows(res_a, res_b2, res_c)

res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])
res <- res %>% drop_na(p)
res <- res %>% arrange(SNP)

res <- res %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
res <- res %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

res <- res %>% group_by(exposure) %>% 
  filter(SNP == first(SNP))

res <- res %>% arrange(p)

res_lung <- res %>% mutate(tissue = "Lung")

res_a <- readRDS("A2_WBC_eQTL_EUR.rds") %>% mutate(outcome = "A2")
res_b2 <- readRDS("B2_WBC_eQTL_EUR.rds") %>% mutate(outcome = "B2")
res_c <- readRDS("C2_WBC_eQTL_EUR.rds") %>% mutate(outcome = "C2")

res <- bind_rows(res_a, res_b2, res_c)

res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])
res <- res %>% drop_na(p)
res <- res %>% arrange(SNP)
res <- res %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
res <- res %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

res <- res %>% group_by(exposure) %>% 
  filter(SNP == first(SNP))


res <- res %>% arrange(p)

res_wbc <- res %>% mutate(tissue = "WBC")

res <- bind_rows(res_lung, res_wbc)


tmp1 <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable1.sQTLMR.lung.xlsx")
tmp1 <- tmp1 %>% filter(p < 0.05/27230)
tmp2 <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable2.sQTLMR.wbc.xlsx")
tmp2 <- tmp2 %>% filter(p < 0.05/27230)

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(res$gene))
colnames(gene) <- "ensembleID"
G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

res <- res %>% left_join(G_list, by=c("gene"="ensembl_gene_id"))

res %>% drop_na(p) %>% filter(gene %in% c(tmp1$gene, tmp2$gene)) %>% filter(p < 0.05)

res <- res %>% mutate(OR = exp(b),
                      LL = exp(b + qnorm(0.025)*se),
                      UL = exp(b + qnorm(0.975)*se),
                      outcome = case_when(outcome == "A2" ~ "Critical illness",
                                          outcome == "B2" ~ "Hospitalization",
                                          outcome == "C2" ~ "Reported infection")
                      ) %>% 
  dplyr::select(gene, hgnc_symbol, outcome, Method, SNPlist, OR, LL, UL, p, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, tissue)
res %>% filter(tissue == "Lung") %>% filter(gene %in% tmp1$gene) 
  write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable6.eQTL_Lung.xlsx")

res %>% filter(tissue == "WBC") %>% filter(gene %in% tmp2$gene) %>% 
  write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable7.eQTL_WBC.xlsx")


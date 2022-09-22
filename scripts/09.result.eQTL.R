setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

#eQTL
res_a <- readRDS("A2_Lung_eQTL.rds") %>% mutate(outcome = "A2")
res_b2 <- readRDS("B2_Lung_eQTL.rds") %>% mutate(outcome = "B2")
res_c <- readRDS("C2_Lung_eQTL.rds") %>% mutate(outcome = "C2")

res <- bind_rows(res_a, res_b2, res_c)

res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])
res <- res %>% drop_na(p)
res <- res %>% arrange(SNP)
res <- res %>% group_by(exposure) %>% 
  filter(SNP == first(SNP))
res <- res %>% mutate(SNP = case_when(grepl("All", SNP) ~ "Inverse variance weighted",
                                      TRUE ~ SNP))

res <- res %>% mutate(p.adj = ifelse(p.adj > 1, 1, p.adj))
res <- res %>% rename(Method = SNP)

res <- res %>% arrange(p)

res_lung <- res %>% mutate(tissue = "Lung")

res_a <- readRDS("A2_WBC_eQTL.rds") %>% mutate(outcome = "A2")
res_b2 <- readRDS("B2_WBC_eQTL.rds") %>% mutate(outcome = "B2")
res_c <- readRDS("C2_WBC_eQTL.rds") %>% mutate(outcome = "C2")

res <- bind_rows(res_a, res_b2, res_c)

res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])

res <- res %>% drop_na(p)
res_a <- readRDS("A2_Lung_eQTL.rds") %>% mutate(outcome = "A2")
res_b2 <- readRDS("B2_Lung_eQTL.rds") %>% mutate(outcome = "B2")
res_c <- readRDS("C2_Lung_eQTL.rds") %>% mutate(outcome = "C2")

res <- bind_rows(res_a, res_b2, res_c)

res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])
res <- res %>% drop_na(p)
res <- res %>% arrange(SNP)
res <- res %>% group_by(exposure) %>% 
  filter(SNP == first(SNP))
res <- res %>% mutate(SNP = case_when(grepl("All", SNP) ~ "Inverse variance weighted",
                                      TRUE ~ SNP))

res <- res %>% mutate(p.adj = ifelse(p.adj > 1, 1, p.adj))
res <- res %>% rename(Method = SNP)
res <- res %>% dplyr::select(gene, exposure, outcome, Method, OR, LL, UL, p, p.adj)

res <- res %>% arrange(p)

res_wbc <- res %>% mutate(tissue = "WBC")

res <- bind_rows(res_lung, res_wbc)


tmp1 <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable1.sQTLMR.lung.xlsx")
tmp1 <- tmp1 %>% filter(p.adj < 0.05)
tmp2 <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable2.sQTLMR.wbc.xlsx")
tmp2 <- tmp2 %>% filter(p.adj < 0.05)

res %>% drop_na(p) %>% filter(gene %in% c(tmp1$gene, tmp2$gene)) %>% filter(p < 0.05)
res %>% drop_na(p) %>% filter(gene %in% c(tmp1$gene, tmp2$gene)) %>% filter(p < 0.05) %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/expression.xlsx")

#res <- res %>% arrange(SNP)
#res <- res %>% group_by(gene) %>% filter(SNP == first(SNP)) %>% ungroup()
res1 <- res %>% dplyr::filter(p < 5e-8)

res <- bind_rows(res_a, res_b2, res_c)
res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])

res <- res %>% dplyr::filter(exposure %in% res1$exposure)

res <- res %>% drop_na(p)

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(res$gene))
colnames(gene) <- "ensembleID"


for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i],mart = ensembl)
}

mrresult <- gene %>% right_join(res, by=c("ensembleID"="gene"))

mrresult %>% write.xlsx("Lung_eQTL.xlsx")

#sQTL
res_a <- readRDS("A2_WBC_eQTL.rds") %>% mutate(outcome = "A2")
res_b2 <- readRDS("B2_WBC_eQTL.rds") %>% mutate(outcome = "B2")
res_c <- readRDS("C2_WBC_eQTL.rds") %>% mutate(outcome = "C2")

res <- bind_rows(res_a, res_b2, res_c)

res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])

res <- res %>% drop_na(p)
#res <- res %>% arrange(SNP)
#res <- res %>% group_by(gene) %>% filter(SNP == first(SNP)) %>% ungroup()
res1 <- res %>% dplyr::filter(p < 5e-8)

res <- bind_rows(res_a, res_b2, res_c)
res <- res %>% dplyr::filter(exposure %in% res1$exposure)

res <- res %>% drop_na(p)
res <- res %>% mutate(gene = str_split(exposure, pattern = "\\.", simplify = T)[,1])


gene <- data.frame(unique(res$gene))
colnames(gene) <- "ensembleID"


for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i],mart = ensembl)
}

mrresult <- gene %>% right_join(res, by=c("ensembleID"="gene"))

mrresult %>% write.xlsx("WBC_eQTL.xlsx")

#eQTL
res_WBC <- readRDS("exposure_eQTL_WBC.rds")
res_Lung <- readRDS("exposure_eQTL_Lung.rds")

mrres <- read.xlsx("WBC_eQTL.xlsx")
res_WBC <- res_WBC %>% dplyr::filter(exposure %in% mrres$exposure)
mrres <- read.xlsx("Lung_eQTL.xlsx")
res_Lung <- res_Lung %>% dplyr::filter(exposure %in% mrres$exposure)

res <- bind_rows(res_WBC, res_Lung) %>% dplyr::select(exposure, chr.exposure, pos.exposure)
res <- res[!duplicated(res$exposure),]
res <- unique(res)
res <- res %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
res <- res %>% dplyr::select(gene, chr.exposure, pos.exposure)
colnames(res) <- c("exposure", "CHR", "POS")
write.table(res, file="eQTL.map", quote=F, col.names = T, row.names = F, sep="\t")

setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

#sQTL
res_a <- readRDS("A2_Lung_sQTL.rds") %>% mutate(outcome = "A2")
res_b2 <- readRDS("B2_Lung_sQTL.rds") %>% mutate(outcome = "B2")
res_c <- readRDS("C2_Lung_sQTL.rds") %>% mutate(outcome = "C2")

res <- bind_rows(res_a, res_b2, res_c)

res <- res %>% mutate(gene = str_split(exposure, pattern = ":", simplify = T)[,5],
                      gene = str_split(gene, pattern = "\\.", simplify = T)[,1])

tmp <- res %>% filter(grepl("chr", SNP))
tmp <- tmp %>%  mutate(CHR = str_split(SNP, pattern = ":",simplify=TRUE)[,1], 
                       POS = as.numeric(str_split(SNP, pattern = ":",simplify=TRUE)[,2]))
tmp <- tmp %>% filter(CHR == "chr6" & POS >= 28510120 & POS <= 33480577)

res <- res %>% filter(!(exposure %in% tmp$exposure))
# chr6:28,510,120-33,480,577

res %>% filter(outcome == "A2") %>% dplyr::select(exposure) %>% unique() %>% dim()#4803
res %>% filter(outcome == "B2") %>% dplyr::select(exposure) %>% unique() %>% dim()#4847
res %>% filter(outcome == "C2") %>% dplyr::select(exposure) %>% unique() %>% dim()#4987
number_test_lung <- 4803+4847+4987
length(unique(res$exposure))#4987
length(unique(res$gene))#3805
res %>% filter(grepl("chr", SNP)) %>% dplyr::select("SNP") %>% unique() %>% dim()#4947
number_test_wbc <- 9373
res <- res %>% mutate(p.adj = ifelse(p*(number_test_wbc + number_test_lung) > 1, 1, p*(number_test_wbc + number_test_lung)))
res %>% filter(p.adj < 0.05) %>% filter(outcome == "A2") %>% dplyr::select(gene) %>% unique() %>% dim()#13
res %>% filter(p.adj < 0.05) %>% filter(outcome == "B2") %>% dplyr::select(gene) %>% unique() %>% dim()#12
res %>% filter(p.adj < 0.05) %>% filter(outcome == "C2") %>% dplyr::select(gene) %>% unique() %>% dim()#11

res <- res %>% mutate(OR = exp(as.numeric(b)),
                      LL = exp(as.numeric(b) - qnorm(0.975)*as.numeric(se)),
                      UL = exp(as.numeric(b) + qnorm(0.975)*as.numeric(se)))

res <- res %>% mutate(outcome = case_when(outcome == "A2" ~ "critical illness",
                                          outcome == "B2" ~ "hospitalization",
                                          TRUE ~ "reported SARS-CoV-2 infection"))
write.xlsx(res, file="All_Lung_sQTL.xlsx")

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

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(res$gene))
colnames(gene) <- "ensembleID"

for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i],mart = ensembl)
}

mrresult <- gene %>% right_join(res, by=c("ensembleID"="gene"))


mrresult %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable1.sQTLMR.lung_rev.xlsx")

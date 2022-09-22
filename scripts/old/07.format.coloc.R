setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

sQTL_WBC <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Whole_Blood.v8.independent_sqtls.txt.gz")
sQTL_WBC <- sQTL_WBC %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
sQTL_WBC <- sQTL_WBC %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                N = round(ma_count/maf*(1/2)))

lung <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable1.sQTLMR.lung.xlsx")
lung <- lung %>% filter(p.adj < 0.05)

wbc <-read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable2.sQTLMR.wbc.xlsx")
wbc <- wbc %>% filter(p.adj < 0.05)

mrresult <- read.xlsx("All_WBC_sQTL.xlsx")
mrresult <- mrresult %>% filter(exposure %in% wbc$exposure)#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult <- mrresult %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

a2 <- fread("A2.coloc.wbc.tsv")
a2 <- a2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
a2 <- a2 %>% mutate(outcome = "critical illness")
b2 <- fread("B2.coloc.wbc.tsv")
b2 <- b2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
b2 <- b2 %>% mutate(outcome = "hospitalization")
c2 <- fread("C2.coloc.wbc.tsv")
c2 <- c2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
c2 <- c2 %>% mutate(outcome = "reported SARS-CoV-2 infection")

colocresult <- bind_rows(a2, b2, c2)

mrresult <- mrresult %>% left_join(colocresult, by=c("junction"="junction","outcome"="outcome"))

mrresult <- mrresult %>% drop_na(p)
mrresult <- mrresult %>% arrange(SNP)
mrresult <- mrresult %>% group_by(exposure, outcome) %>%
  filter(SNP == first(SNP))
  
library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(mrresult$gene))
colnames(gene) <- "ensembleID"

for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i],mart = ensembl)
}

mrresult <- gene %>% right_join(mrresult, by=c("ensembleID"="gene"))

write.xlsx(mrresult, file="WBC.coloc_sQTL.xlsx")

mrresult <- mrresult %>% dplyr::select(ensembleID, gene_id, phenotype_id,
                                PP.H0.abf, PP.H1.abf, PP.H2.abf,
                                PP.H3.abf, PP.H4.abf, exposure, outcome,
                                SNP, OR, LL, UL, p, p.adj)

mrresult %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable4.coloc.wbc.xlsx")

sQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Lung.v8.independent_sqtls.txt.gz")
sQTL_Lung <- sQTL_Lung %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
sQTL_Lung <- sQTL_Lung %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                N = round(ma_count/maf*(1/2)))

#sQTL_Lung %>% filter(grepl("chr13:112875941:112880546:clu_3196:ENSG00000068650.18", phenotype_id))
mrresult <- read.xlsx("All_Lung_sQTL.xlsx")
mrresult <- mrresult %>% filter(exposure %in% lung$exposure)#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult <- mrresult %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                  str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                  str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                  str_split(exposure, pattern = ":", simplify = T)[,5]))

a2 <- fread("A2.coloc.lung.tsv")
a2 <- a2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
a2 <- a2 %>% mutate(outcome = "critical illness")
b2 <- fread("B2.coloc.lung.tsv")
b2 <- b2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
b2 <- b2 %>% mutate(outcome = "hospitalization")
c2 <- fread("C2.coloc.lung.tsv")
c2 <- c2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
c2 <- c2 %>% mutate(outcome = "reported SARS-CoV-2 infection")

colocresult <- bind_rows(a2, b2, c2)

mrresult <- mrresult %>% left_join(colocresult, by=c("junction"="junction","outcome"="outcome"))

mrresult <- mrresult %>% drop_na(p)
mrresult <- mrresult %>% arrange(SNP)
mrresult <- mrresult %>% group_by(exposure, outcome) %>%
  filter(SNP == first(SNP))

gene <- data.frame(unique(mrresult$gene))
colnames(gene) <- "ensembleID"

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(mrresult$gene))
colnames(gene) <- "ensembleID"


for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i],mart = ensembl)
}

mrresult <- gene %>% right_join(mrresult, by=c("ensembleID"="gene"))

write.xlsx(mrresult, file="Lung.coloc_sQTL.xlsx")

mrresult <- mrresult %>% dplyr::select(ensembleID, gene_id, phenotype_id,
                                PP.H0.abf, PP.H1.abf, PP.H2.abf,
                                PP.H3.abf, PP.H4.abf, exposure, outcome,
                                SNP, OR, LL, UL, p, p.adj)

mrresult %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable3.coloc.lung.xlsx")


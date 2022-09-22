setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

coloc_result <- read.xlsx("Lung.coloc.xlsx")

coloc_result <- coloc_result %>% filter(PP.H4.abf > PP.H0.abf + PP.H1.abf + PP.H2.abf + PP.H3.abf)
coloc_result <- coloc_result %>%  mutate(CHR = str_split(SNP, pattern = ":",simplify=TRUE)[,1], 
                                         POS = as.numeric(str_split(SNP, pattern = ":",simplify=TRUE)[,2]))

eQTL_WBC <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_Analysis_v8_eQTL_independent/Whole_Blood.v8.independent_eqtls.txt.gz")
eQTL_WBC <- eQTL_WBC %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
eQTL_WBC <- eQTL_WBC %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                N = round(minor_allele_count/maf*(1/2)))
eQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_Analysis_v8_eQTL_independent/Lung.v8.independent_eqtls.txt.gz")
eQTL_Lung <- eQTL_Lung %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                  POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                  EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                  NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
eQTL_Lung <- eQTL_Lung %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                  N = round(minor_allele_count/maf*(1/2)))

#ATP11A
a_wbc_mr <- readRDS("A2_WBC_eQTL.rds")
a_wbc_mr <- a_wbc_mr %>% filter(grepl(coloc_result$ensembleID[12], exposure))
a_wbc_mr <- a_wbc_mr %>% drop_na(p)
a_wbc_coloc <- fread("A2.coloc.wbc.eQTL.tsv")
a_wbc_mr <- a_wbc_mr %>% inner_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
a_wbc_mr <- a_wbc_mr %>% mutate(outcome = "A2") %>% mutate(tissue = "WBC")

b_wbc_mr <- readRDS("B2_WBC_eQTL.rds")
b_wbc_mr <- b_wbc_mr %>% filter(grepl(coloc_result$ensembleID[12], exposure))
b_wbc_mr <- b_wbc_mr %>% drop_na(p)
b_wbc_coloc <- fread("B2.coloc.wbc.eQTL.tsv")
b_wbc_mr <- b_wbc_mr %>% inner_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
b_wbc_mr <- b_wbc_mr %>% mutate(outcome = "B2") %>% mutate(tissue = "WBC")

c_wbc_mr <- readRDS("C2_WBC_eQTL.rds")
c_wbc_mr <- c_wbc_mr %>% filter(grepl(coloc_result$ensembleID[12], exposure))
c_wbc_mr <- c_wbc_mr %>% drop_na(p)
c_wbc_coloc <- fread("C2.coloc.wbc.eQTL.tsv")
c_wbc_mr <- c_wbc_mr %>% inner_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
c_wbc_mr <- c_wbc_mr %>% mutate(outcome = "C2") %>% mutate(tissue = "WBC")

a_lung_coloc <- fread("A2.coloc.Lung.eQTL.tsv")
a_lung_coloc <- a_lung_coloc %>% filter(grepl(coloc_result$ensembleID[12], phenotype_id))
a_lung_coloc <- a_lung_coloc %>% mutate(outcome = "A2") %>% mutate(tissue = "lung")

b_lung_coloc <- fread("B2.coloc.Lung.eQTL.tsv")
b_lung_coloc <- b_lung_coloc %>% filter(grepl(coloc_result$ensembleID[12], phenotype_id))
b_lung_coloc <- b_lung_coloc %>% mutate(outcome = "B2") %>% mutate(tissue = "lung")

c_lung_coloc <- fread("C2.coloc.Lung.eQTL.tsv")
c_lung_coloc <- c_lung_coloc %>% filter(grepl(coloc_result$ensembleID[12], phenotype_id))
c_lung_coloc <- c_lung_coloc %>% mutate(outcome = "C2") %>% mutate(tissue = "lung")


ATP11A_result <- bind_rows(a_wbc_mr, b_wbc_mr, c_wbc_mr,  a_lung_coloc, 
                           b_lung_coloc, c_lung_coloc)

write.xlsx(c_mr, "ATP11A_expression.xlsx")

#DPP9
a_wbc_mr <- readRDS("A2_WBC_eQTL.rds")
a_wbc_mr <- a_wbc_mr %>% filter(grepl(coloc_result$ensembleID[4], exposure))
a_wbc_mr <- a_wbc_mr %>% drop_na(p)
a_wbc_coloc <- fread("A2.coloc.wbc.eQTL.tsv")
a_wbc_mr <- a_wbc_mr %>% inner_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
a_wbc_mr <- a_wbc_mr %>% mutate(outcome = "A2") %>% mutate(tissue = "WBC")

b_wbc_mr <- readRDS("B2_WBC_eQTL.rds")
b_wbc_mr <- b_wbc_mr %>% filter(grepl(coloc_result$ensembleID[4], exposure))
b_wbc_mr <- b_wbc_mr %>% drop_na(p)
b_wbc_coloc <- fread("B2.coloc.wbc.eQTL.tsv")
b_wbc_mr <- b_wbc_mr %>% inner_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
b_wbc_mr <- b_wbc_mr %>% mutate(outcome = "B2") %>% mutate(tissue = "WBC")

c_wbc_mr <- readRDS("C2_WBC_eQTL.rds")
c_wbc_mr <- c_wbc_mr %>% filter(grepl(coloc_result$ensembleID[4], exposure))
c_wbc_mr <- c_wbc_mr %>% drop_na(p)
c_wbc_coloc <- fread("C2.coloc.wbc.eQTL.tsv")
c_wbc_mr <- c_wbc_mr %>% inner_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
c_wbc_mr <- c_wbc_mr %>% mutate(outcome = "C2") %>% mutate(tissue = "WBC")

a_lung_coloc <- fread("A2.coloc.Lung.eQTL.tsv")
a_lung_coloc <- a_lung_coloc %>% filter(grepl(coloc_result$ensembleID[4], phenotype_id))
a_lung_coloc <- a_lung_coloc %>% mutate(outcome = "A2") %>% mutate(tissue = "lung")

b_lung_coloc <- fread("B2.coloc.Lung.eQTL.tsv")
b_lung_coloc <- b_lung_coloc %>% filter(grepl(coloc_result$ensembleID[4], phenotype_id))
b_lung_coloc <- b_lung_coloc %>% mutate(outcome = "B2") %>% mutate(tissue = "lung")

c_lung_coloc <- fread("C2.coloc.Lung.eQTL.tsv")
c_lung_coloc <- c_lung_coloc %>% filter(grepl(coloc_result$ensembleID[4], phenotype_id))
c_lung_coloc <- c_lung_coloc %>% mutate(outcome = "C2") %>% mutate(tissue = "lung")


DPP9_result <- bind_rows(a_wbc_mr, b_wbc_mr, c_wbc_mr,  a_lung_coloc, 
                           b_lung_coloc, c_lung_coloc)

##NPNT
a_wbc_coloc <- fread("A2.coloc.wbc.eQTL.tsv")
a_wbc_coloc <- a_wbc_coloc %>% filter(grepl(coloc_result$ensembleID[7], phenotype_id))
a_wbc_coloc <- a_wbc_coloc %>% mutate(outcome = "A2") %>% mutate(tissue = "wbc")

b_wbc_coloc <- fread("B2.coloc.wbc.eQTL.tsv")
b_wbc_coloc <- b_wbc_coloc %>% filter(grepl(coloc_result$ensembleID[7], phenotype_id))
b_wbc_coloc <- b_wbc_coloc %>% mutate(outcome = "B2") %>% mutate(tissue = "wbc")

c_wbc_coloc <- fread("C2.coloc.wbc.eQTL.tsv")
c_wbc_coloc <- c_wbc_coloc %>% filter(grepl(coloc_result$ensembleID[7], phenotype_id))
c_wbc_coloc <- c_wbc_coloc %>% mutate(outcome = "C2") %>% mutate(tissue = "wbc")

a_lung_mr <- readRDS("A2_Lung_eQTL.rds")
a_lung_mr <- a_lung_mr %>% filter(grepl(coloc_result$ensembleID[7], exposure))
a_lung_mr <- a_lung_mr %>% filter(grepl("All - Inverse", SNP))
a_lung_coloc <- fread("A2.coloc.Lung.eQTL.tsv")
a_lung_mr <- a_lung_mr %>% inner_join(a_lung_coloc, by=c("exposure"="phenotype_id"))
a_lung_mr <- a_lung_mr %>% mutate(outcome = "A2") %>% mutate(tissue = "lung")

b_lung_mr <- readRDS("B2_Lung_eQTL.rds")
b_lung_mr <- b_lung_mr %>% filter(grepl(coloc_result$ensembleID[7], exposure))
b_lung_mr <- b_lung_mr %>% filter(grepl("All - Inverse", SNP))
b_lung_coloc <- fread("B2.coloc.Lung.eQTL.tsv")
b_lung_mr <- b_lung_mr %>% inner_join(a_lung_coloc, by=c("exposure"="phenotype_id"))
b_lung_mr <- b_lung_mr %>% mutate(outcome = "B2") %>% mutate(tissue = "lung")

c_lung_mr <- readRDS("C2_Lung_eQTL.rds")
c_lung_mr <- c_lung_mr %>% filter(grepl(coloc_result$ensembleID[7], exposure))
c_lung_mr <- c_lung_mr %>% filter(grepl("All - Inverse", SNP))
c_lung_coloc <- fread("C2.coloc.Lung.eQTL.tsv")
c_lung_mr <- c_lung_mr %>% inner_join(a_lung_coloc, by=c("exposure"="phenotype_id"))
c_lung_mr <- c_lung_mr %>% mutate(outcome = "C2") %>% mutate(tissue = "lung")

NPNT_result <- bind_rows(a_lung_mr, b_lung_mr, c_lung_mr,  a_wbc_coloc, 
                         b_wbc_coloc, c_wbc_coloc)

##SFTPA2
a_lung_mr <- readRDS("A2_Lung_eQTL.rds")
a_lung_mr <- a_lung_mr %>% filter(grepl(coloc_result$ensembleID[9], exposure))
a_lung_mr <- a_lung_mr %>% drop_na(p)
a_lung_coloc <- fread("A2.coloc.Lung.eQTL.tsv")
a_lung_mr <- a_lung_mr %>% inner_join(a_lung_coloc, by=c("exposure"="phenotype_id"))
a_lung_mr <- a_lung_mr %>% mutate(outcome = "A2") %>% mutate(tissue = "lung")

b_lung_mr <- readRDS("B2_Lung_eQTL.rds")
b_lung_mr <- b_lung_mr %>% filter(grepl(coloc_result$ensembleID[9], exposure))
b_lung_mr <- b_lung_mr %>% drop_na(p)
b_lung_coloc <- fread("B2.coloc.Lung.eQTL.tsv")
b_lung_mr <- b_lung_mr %>% inner_join(b_lung_coloc, by=c("exposure"="phenotype_id"))
b_lung_mr <- b_lung_mr %>% mutate(outcome = "B2") %>% mutate(tissue = "lung")

c_lung_mr <- readRDS("C2_Lung_eQTL.rds")
c_lung_mr <- c_lung_mr %>% filter(grepl(coloc_result$ensembleID[9], exposure))
c_lung_mr <- c_lung_mr %>% drop_na(p)
c_lung_coloc <- fread("C2.coloc.Lung.eQTL.tsv")
c_lung_mr <- c_lung_mr %>% inner_join(c_lung_coloc, by=c("exposure"="phenotype_id"))
c_lung_mr <- c_lung_mr %>% mutate(outcome = "C2") %>% mutate(tissue = "lung")

SFTPA2_result <- bind_rows(a_lung_mr, b_lung_mr, c_lung_mr)

##SFTPA1
a_lung_mr <- readRDS("A2_Lung_eQTL.rds")
a_lung_mr <- a_lung_mr %>% filter(grepl("ENSG00000185303", exposure))
a_lung_mr <- a_lung_mr %>% drop_na(p)
a_lung_coloc <- fread("A2.coloc.Lung.eQTL.tsv")
a_lung_mr <- a_lung_mr %>% inner_join(a_lung_coloc, by=c("exposure"="phenotype_id"))
a_lung_mr <- a_lung_mr %>% mutate(outcome = "A2") %>% mutate(tissue = "lung")

b_lung_mr <- readRDS("B2_Lung_eQTL.rds")
b_lung_mr <- b_lung_mr %>% filter(grepl("ENSG00000185303", exposure))
b_lung_mr <- b_lung_mr %>% drop_na(p)
b_lung_coloc <- fread("B2.coloc.Lung.eQTL.tsv")
b_lung_mr <- b_lung_mr %>% inner_join(b_lung_coloc, by=c("exposure"="phenotype_id"))
b_lung_mr <- b_lung_mr %>% mutate(outcome = "B2") %>% mutate(tissue = "lung")

c_lung_mr <- readRDS("C2_Lung_eQTL.rds")
c_lung_mr <- c_lung_mr %>% filter(grepl("ENSG00000185303", exposure))
c_lung_mr <- c_lung_mr %>% drop_na(p)
c_lung_coloc <- fread("C2.coloc.Lung.eQTL.tsv")
c_lung_mr <- c_lung_mr %>% inner_join(c_lung_coloc, by=c("exposure"="phenotype_id"))
c_lung_mr <- c_lung_mr %>% mutate(outcome = "C2") %>% mutate(tissue = "lung")

SFTPA2_result <- bind_rows(a_lung_mr, b_lung_mr, c_lung_mr)




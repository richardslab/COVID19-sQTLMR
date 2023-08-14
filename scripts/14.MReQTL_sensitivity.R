setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")
library(tidyverse)
library(data.table)
library(openxlsx)

ALL_Lung_sQTL_EUR <- readRDS("EUR_Lung_sQTL_EUR_coloc.rds")
ALL_WBC_sQTL_EUR <- readRDS("EUR_WBC_sQTL_EUR_coloc.rds")

# tmp_lung <- ALL_Lung_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8)
# tmp_wbc <- ALL_WBC_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8)
# 
# tmp <- bind_rows(tmp_lung, tmp_wbc)

eQTL_list <- c("THBS3", "EFNA1", "ADAM15", "GBA", "MTX1", "YY1AP1", "DAP3", 
               "TRIM46", "FAM189B", "HCN3", "RIT1", "DPM3", "KRTCAP2", "SYT11", 
               "SEMA4A", "PKLR", "FDPS", "SCAMP3", "MUC1", "GBAP1", 
               "GLMP", "SMG5", "PAQR6", "TMEME79", "SLC25A44", "CCT3",
               "OAS3", "RPH3A", "OAS2", "OAS1",
               "TNFAIP8L1", "PMF1", "NPNT","ATP11A","DPP9","ABO","GBGT1")
eQTL_list <- unique(eQTL_list)
library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(eQTL_list)
colnames(gene) <- "symbol"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'hgnc_symbol',values = gene$symbol,mart = ensembl)


A2_Lung_eQTL_EUR <- readRDS("A2_Lung_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% drop_na(p) %>% arrange(SNP)
A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% group_by(exposure) %>%
  filter(SNP == first(SNP))
A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% filter(gene %in% G_list$ensembl_gene_id)

A2_Lung_eQTL_coloc <- fread("GTEx/eQTL/A2_eQTL_GTEx_Lung_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
A2_Lung_eQTL_coloc <- A2_Lung_eQTL_coloc %>% filter(gene %in% G_list$ensembl_gene_id)
A2_Lung_eQTL_coloc <- A2_Lung_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
A2_Lung <- A2_Lung_eQTL_EUR %>% full_join(A2_Lung_eQTL_coloc, by=c("gene"="gene"))

A2_Lung <- A2_Lung %>% mutate(outcome = "Critical illness")

###########
B2_Lung_eQTL_EUR <- readRDS("B2_Lung_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
B2_Lung_eQTL_EUR <- B2_Lung_eQTL_EUR %>% drop_na(p) %>% arrange(SNP)
B2_Lung_eQTL_EUR <- B2_Lung_eQTL_EUR %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
B2_Lung_eQTL_EUR <- B2_Lung_eQTL_EUR %>% group_by(exposure) %>%
  filter(SNP == first(SNP))
B2_Lung_eQTL_EUR <- B2_Lung_eQTL_EUR %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

B2_Lung_eQTL_EUR <- B2_Lung_eQTL_EUR %>% filter(gene %in% G_list$ensembl_gene_id)


B2_Lung_eQTL_coloc <- fread("GTEx/eQTL/B2_eQTL_GTEx_Lung_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
B2_Lung_eQTL_coloc <- B2_Lung_eQTL_coloc %>% filter(gene %in% G_list$ensembl_gene_id)
B2_Lung_eQTL_coloc <- B2_Lung_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
B2_Lung <- B2_Lung_eQTL_EUR %>% full_join(B2_Lung_eQTL_coloc, by=c("gene"="gene"))

B2_Lung <- B2_Lung %>% mutate(outcome = "Hospitalization")

###############
C2_Lung_eQTL_EUR <- readRDS("C2_Lung_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
C2_Lung_eQTL_EUR <- C2_Lung_eQTL_EUR %>% drop_na(p) %>% arrange(SNP)
C2_Lung_eQTL_EUR <- C2_Lung_eQTL_EUR %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
C2_Lung_eQTL_EUR <- C2_Lung_eQTL_EUR %>% group_by(exposure) %>%
  filter(SNP == first(SNP))
C2_Lung_eQTL_EUR <- C2_Lung_eQTL_EUR %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

C2_Lung_eQTL_EUR <- C2_Lung_eQTL_EUR %>% filter(gene %in% G_list$ensembl_gene_id)

C2_Lung_eQTL_coloc <- fread("GTEx/eQTL/C2_eQTL_GTEx_Lung_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
C2_Lung_eQTL_coloc <- C2_Lung_eQTL_coloc %>% filter(gene %in% G_list$ensembl_gene_id)
C2_Lung_eQTL_coloc <- C2_Lung_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
C2_Lung <- C2_Lung_eQTL_EUR %>% full_join(C2_Lung_eQTL_coloc, by=c("gene"="gene"))

C2_Lung <- C2_Lung %>% mutate(outcome = "Reported infection")

##################

Lung <- bind_rows(A2_Lung, B2_Lung, C2_Lung)

Lung <- Lung %>% mutate(OR = exp(b),
                        LL = exp(b + qnorm(0.025)*se),
                        UL = exp(b + qnorm(0.975)*se)) %>% 
  dplyr::select(gene, outcome, Method, SNPlist, OR, LL, UL, p, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

Lung <- Lung %>% left_join(G_list, by=c("gene"="ensembl_gene_id"))

Lung %>% filter(PP.H4.abf > 0.8)
Lung %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable7.eQTL_Lung.xlsx")


################WBC
A2_WBC_eQTL_EUR <- readRDS("A2_WBC_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
A2_WBC_eQTL_EUR <- A2_WBC_eQTL_EUR %>% drop_na(p) %>% arrange(SNP)
A2_WBC_eQTL_EUR <- A2_WBC_eQTL_EUR %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
A2_WBC_eQTL_EUR <- A2_WBC_eQTL_EUR %>% group_by(exposure) %>%
  filter(SNP == first(SNP))
A2_WBC_eQTL_EUR <- A2_WBC_eQTL_EUR %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

A2_WBC_eQTL_EUR <- A2_WBC_eQTL_EUR %>% filter(gene %in% G_list$ensembl_gene_id)

A2_WBC_eQTL_coloc <- fread("GTEx/eQTL/A2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
A2_WBC_eQTL_coloc <- A2_WBC_eQTL_coloc %>% filter(gene %in% G_list$ensembl_gene_id)
A2_WBC_eQTL_coloc <- A2_WBC_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
A2_WBC <- A2_WBC_eQTL_EUR %>% full_join(A2_WBC_eQTL_coloc, by=c("gene"="gene"))

A2_WBC <- A2_WBC %>% mutate(outcome = "Critical illness")

###########
B2_WBC_eQTL_EUR <- readRDS("B2_WBC_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
B2_WBC_eQTL_EUR <- B2_WBC_eQTL_EUR %>% drop_na(p) %>% arrange(SNP)
B2_WBC_eQTL_EUR <- B2_WBC_eQTL_EUR %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
B2_WBC_eQTL_EUR <- B2_WBC_eQTL_EUR %>% group_by(exposure) %>%
  filter(SNP == first(SNP))
B2_WBC_eQTL_EUR <- B2_WBC_eQTL_EUR %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

B2_WBC_eQTL_EUR <- B2_WBC_eQTL_EUR %>% filter(gene %in% G_list$ensembl_gene_id)

B2_WBC_eQTL_coloc <- fread("GTEx/eQTL/B2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
B2_WBC_eQTL_coloc <- B2_WBC_eQTL_coloc %>% filter(gene %in% G_list$ensembl_gene_id)
B2_WBC_eQTL_coloc <- B2_WBC_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
B2_WBC <- B2_WBC_eQTL_EUR %>% full_join(B2_WBC_eQTL_coloc, by=c("gene"="gene"))

B2_WBC <- B2_WBC %>% mutate(outcome = "Hospitalization")

###############
C2_WBC_eQTL_EUR <- readRDS("C2_WBC_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
C2_WBC_eQTL_EUR <- C2_WBC_eQTL_EUR %>% drop_na(p) %>% arrange(SNP)
C2_WBC_eQTL_EUR <- C2_WBC_eQTL_EUR %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
C2_WBC_eQTL_EUR <- C2_WBC_eQTL_EUR %>% group_by(exposure) %>%
  filter(SNP == first(SNP))
C2_WBC_eQTL_EUR <- C2_WBC_eQTL_EUR %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

C2_WBC_eQTL_EUR <- C2_WBC_eQTL_EUR %>% filter(gene %in% G_list$ensembl_gene_id)

C2_WBC_eQTL_coloc <- fread("GTEx/eQTL/C2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
C2_WBC_eQTL_coloc <- C2_WBC_eQTL_coloc %>% filter(gene %in% G_list$ensembl_gene_id)
C2_WBC_eQTL_coloc <- C2_WBC_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
C2_WBC <- C2_WBC_eQTL_EUR %>% full_join(C2_WBC_eQTL_coloc, by=c("gene"="gene"))

C2_WBC <- C2_WBC %>% mutate(outcome = "Reported infection")

##################

WBC <- bind_rows(A2_WBC, B2_WBC, C2_WBC)


WBC <- WBC %>% mutate(OR = exp(b),
                      LL = exp(b + qnorm(0.025)*se),
                      UL = exp(b + qnorm(0.975)*se)) %>% 
  dplyr::select(gene, outcome, Method, SNPlist, OR, LL, UL, p, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

WBC <- WBC %>% left_join(G_list, by=c("gene"="ensembl_gene_id"))

WBC %>% filter(PP.H4.abf > 0.8) %>% View()
WBC %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable8.eQTL_WBC.xlsx")


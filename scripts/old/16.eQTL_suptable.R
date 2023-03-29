setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")
library(tidyverse)
library(data.table)
library(openxlsx)

ALL_Lung_sQTL_EUR <- readRDS("EUR_Lung_sQTL_EUR_coloc.rds")
ALL_WBC_sQTL_EUR <- readRDS("EUR_WBC_sQTL_EUR_coloc.rds")

tmp_lung <- ALL_Lung_sQTL_EUR %>% filter(p < 0.05/9292 & PP.H4.abf > 0.8)
tmp_wbc <- ALL_WBC_sQTL_EUR %>% filter(p < 0.05/9292 & PP.H4.abf > 0.8)

tmp <- bind_rows(tmp_lung, tmp_wbc)

A2_Lung_eQTL_ALL <- readRDS("A2_Lung_eQTL_ALL.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
A2_Lung_eQTL_ALL <- A2_Lung_eQTL_ALL %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
A2_Lung_eQTL_ALL <- A2_Lung_eQTL_ALL %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

A2_Lung_eQTL_ALL <- A2_Lung_eQTL_ALL %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)

A2_Lung_eQTL_EUR <- readRDS("A2_Lung_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=",")) %>%
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))

A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
A2_Lung_eQTL_EUR <- A2_Lung_eQTL_EUR %>% dplyr::select(gene, SNP, b, se, p)
A2_Lung_eQTL_coloc <- fread("GTEx/eQTL/A2_eQTL_GTEx_Lung_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
A2_Lung_eQTL_coloc <- A2_Lung_eQTL_coloc %>% filter(gene %in% tmp$ensembleID)
A2_Lung_eQTL_coloc <- A2_Lung_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

length(unique(A2_Lung_eQTL_ALL$gene))
length(unique(A2_Lung_eQTL_EUR$gene))
length(unique(A2_Lung_eQTL_coloc$gene))
A2_Lung <- A2_Lung_eQTL_ALL %>% full_join(A2_Lung_eQTL_EUR, by=c("gene"="gene", "SNP"="SNP"))
A2_Lung <- A2_Lung %>% full_join(A2_Lung_eQTL_coloc, by=c("gene"="gene"))

A2_Lung <- A2_Lung %>% mutate(outcome = "Critical illness")

B2_Lung_eQTL_ALL <- readRDS("B2_Lung_eQTL_ALL.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
B2_Lung_eQTL_ALL <- B2_Lung_eQTL_ALL %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
B2_Lung_eQTL_ALL <- B2_Lung_eQTL_ALL %>% dplyr::select(gene, SNP, b, se, p)
B2_Lung_eQTL_EUR <- readRDS("B2_Lung_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
B2_Lung_eQTL_EUR <- B2_Lung_eQTL_EUR %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
B2_Lung_eQTL_EUR <- B2_Lung_eQTL_EUR %>% dplyr::select(gene, SNP, b, se, p)
B2_Lung_eQTL_coloc <- fread("GTEx/eQTL/B2_eQTL_GTEx_Lung_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
B2_Lung_eQTL_coloc <- B2_Lung_eQTL_coloc %>% filter(gene %in% tmp$ensembleID)
B2_Lung_eQTL_coloc <- B2_Lung_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
length(unique(B2_Lung_eQTL_ALL$gene))
length(unique(B2_Lung_eQTL_EUR$gene))
length(unique(B2_Lung_eQTL_coloc$gene))
B2_Lung <- B2_Lung_eQTL_ALL %>% full_join(B2_Lung_eQTL_EUR, by=c("gene"="gene", "SNP"="SNP"))
B2_Lung <- B2_Lung %>% full_join(B2_Lung_eQTL_coloc, by=c("gene"="gene"))

B2_Lung <- B2_Lung %>% mutate(outcome = "Hospitalization")

C2_Lung_eQTL_ALL <- readRDS("C2_Lung_eQTL_ALL.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
C2_Lung_eQTL_ALL <- C2_Lung_eQTL_ALL %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
C2_Lung_eQTL_ALL <- C2_Lung_eQTL_ALL %>% dplyr::select(gene, SNP, b, se, p)
C2_Lung_eQTL_EUR <- readRDS("C2_Lung_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
C2_Lung_eQTL_EUR <- C2_Lung_eQTL_EUR %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
C2_Lung_eQTL_EUR <- C2_Lung_eQTL_EUR %>% dplyr::select(gene, SNP, b, se, p)
C2_Lung_eQTL_coloc <- fread("GTEx/eQTL/C2_eQTL_GTEx_Lung_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
C2_Lung_eQTL_coloc <- C2_Lung_eQTL_coloc %>% filter(gene %in% tmp$ensembleID)
C2_Lung_eQTL_coloc <- C2_Lung_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
length(unique(C2_Lung_eQTL_ALL$gene))
length(unique(C2_Lung_eQTL_EUR$gene))
length(unique(C2_Lung_eQTL_coloc$gene))
C2_Lung <- C2_Lung_eQTL_ALL %>% full_join(C2_Lung_eQTL_EUR, by=c("gene"="gene", "SNP"="SNP"))
C2_Lung <- C2_Lung %>% full_join(C2_Lung_eQTL_coloc, by=c("gene"="gene"))

C2_Lung <- C2_Lung %>% mutate(outcome = "Reported infection")

Lung <- bind_rows(A2_Lung, B2_Lung, C2_Lung)
Lung <- Lung %>% mutate(OR.ALL = exp(b.x),
                        LL.ALL = exp(b.x + qnorm(0.025)*se.x),
                        UL.ALL = exp(b.x + qnorm(0.975)*se.x),
                        p.ALL = p.x,
                        OR.EUR = exp(b.y),
                        LL.EUR = exp(b.y + qnorm(0.025)*se.y),
                        UL.EUR = exp(b.y + qnorm(0.975)*se.y),
                        p.EUR = p.y)

Lung <- Lung %>% dplyr::select(gene, outcome, SNP, OR.ALL, LL.ALL, UL.ALL, p.ALL, OR.EUR, LL.EUR, UL.EUR,
                               p.EUR, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(Lung$gene))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

Lung <- Lung %>% left_join(G_list, by=c("gene"="ensembl_gene_id"))

Lung <- Lung %>% dplyr::select(gene, hgnc_symbol, outcome, SNP, OR.ALL, LL.ALL,UL.ALL, p.ALL, OR.EUR,
                               LL.EUR,UL.EUR, p.EUR,NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

Lung %>% filter(PP.H4.abf > 0.8)
Lung %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable5.eQTL_Lung.xlsx")


A2_WBC_eQTL_ALL <- readRDS("A2_WBC_eQTL_ALL.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
A2_WBC_eQTL_ALL <- A2_WBC_eQTL_ALL %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
A2_WBC_eQTL_ALL <- A2_WBC_eQTL_ALL %>% dplyr::select(gene, SNP, b, se, p)
A2_WBC_eQTL_EUR <- readRDS("A2_WBC_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
A2_WBC_eQTL_EUR <- A2_WBC_eQTL_EUR %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
A2_WBC_eQTL_EUR <- A2_WBC_eQTL_EUR %>% dplyr::select(gene, SNP, b, se, p)
A2_WBC_eQTL_coloc <- fread("GTEx/eQTL/A2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
A2_WBC_eQTL_coloc <- A2_WBC_eQTL_coloc %>% filter(gene %in% tmp$ensembleID)
A2_WBC_eQTL_coloc <- A2_WBC_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
A2_WBC_eQTL_coloc_eQTLgen <- fread("eQTLgen/A2_eQTL_eQTLgen.tsv") %>% mutate(gene = phenotype_id)
A2_WBC_eQTL_coloc_eQTLgen <- A2_WBC_eQTL_coloc_eQTLgen %>% filter(gene %in% tmp$ensembleID)
A2_WBC_eQTL_coloc_eQTLgen <- A2_WBC_eQTL_coloc_eQTLgen %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

length(unique(A2_WBC_eQTL_ALL$gene))
length(unique(A2_WBC_eQTL_EUR$gene))
length(unique(A2_WBC_eQTL_coloc$gene))
length(unique(A2_WBC_eQTL_coloc_eQTLgen$gene))
A2_WBC <- A2_WBC_eQTL_ALL %>% full_join(A2_WBC_eQTL_EUR, by=c("gene"="gene", "SNP"="SNP"))
A2_WBC <- A2_WBC %>% full_join(A2_WBC_eQTL_coloc, by=c("gene"="gene"))
A2_WBC <- A2_WBC %>% full_join(A2_WBC_eQTL_coloc_eQTLgen, by=c("gene"="gene"))

A2_WBC <- A2_WBC %>% mutate(outcome = "Critical illness")

B2_WBC_eQTL_ALL <- readRDS("B2_WBC_eQTL_ALL.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
B2_WBC_eQTL_ALL <- B2_WBC_eQTL_ALL %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
B2_WBC_eQTL_ALL <- B2_WBC_eQTL_ALL %>% dplyr::select(gene, SNP, b, se, p)
B2_WBC_eQTL_EUR <- readRDS("B2_WBC_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
B2_WBC_eQTL_EUR <- B2_WBC_eQTL_EUR %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
B2_WBC_eQTL_EUR <- B2_WBC_eQTL_EUR %>% dplyr::select(gene, SNP, b, se, p)
B2_WBC_eQTL_coloc <- fread("GTEx/eQTL/B2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
B2_WBC_eQTL_coloc <- B2_WBC_eQTL_coloc %>% filter(gene %in% tmp$ensembleID)
B2_WBC_eQTL_coloc <- B2_WBC_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
B2_WBC_eQTL_coloc_eQTLgen <- fread("eQTLgen/B2_eQTL_eQTLgen.tsv") %>% mutate(gene = phenotype_id)
B2_WBC_eQTL_coloc_eQTLgen <- B2_WBC_eQTL_coloc_eQTLgen %>% filter(gene %in% tmp$ensembleID)
B2_WBC_eQTL_coloc_eQTLgen <- B2_WBC_eQTL_coloc_eQTLgen %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)


length(unique(B2_WBC_eQTL_ALL$gene))
length(unique(B2_WBC_eQTL_EUR$gene))
length(unique(B2_WBC_eQTL_coloc$gene))
length(unique(B2_WBC_eQTL_coloc_eQTLgen$gene))
B2_WBC <- B2_WBC_eQTL_ALL %>% full_join(B2_WBC_eQTL_EUR, by=c("gene"="gene", "SNP"="SNP"))
B2_WBC <- B2_WBC %>% full_join(B2_WBC_eQTL_coloc, by=c("gene"="gene"))
B2_WBC <- B2_WBC %>% full_join(B2_WBC_eQTL_coloc_eQTLgen, by=c("gene"="gene"))
B2_WBC <- B2_WBC %>% mutate(outcome = "Hospitalization")

C2_WBC_eQTL_ALL <- readRDS("C2_WBC_eQTL_ALL.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
C2_WBC_eQTL_ALL <- C2_WBC_eQTL_ALL %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
C2_WBC_eQTL_ALL <- C2_WBC_eQTL_ALL %>% dplyr::select(gene, SNP, b, se, p)
C2_WBC_eQTL_EUR <- readRDS("C2_WBC_eQTL_EUR.rds") %>% mutate(gene = str_split(exposure, pattern="\\.", simplify = T)[,1])
C2_WBC_eQTL_EUR <- C2_WBC_eQTL_EUR %>% filter(gene %in% tmp$ensembleID) %>% drop_na(p)
C2_WBC_eQTL_EUR <- C2_WBC_eQTL_EUR %>% dplyr::select(gene, SNP, b, se, p)
C2_WBC_eQTL_coloc <- fread("GTEx/eQTL/C2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz") %>% mutate(gene = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
C2_WBC_eQTL_coloc <- C2_WBC_eQTL_coloc %>% filter(gene %in% tmp$ensembleID)
C2_WBC_eQTL_coloc <- C2_WBC_eQTL_coloc %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
C2_WBC_eQTL_coloc_eQTLgen <- fread("eQTLgen/C2_eQTL_eQTLgen.tsv") %>% mutate(gene = phenotype_id)
C2_WBC_eQTL_coloc_eQTLgen <- C2_WBC_eQTL_coloc_eQTLgen %>% filter(gene %in% tmp$ensembleID)
C2_WBC_eQTL_coloc_eQTLgen <- C2_WBC_eQTL_coloc_eQTLgen %>% dplyr::select(gene, NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

length(unique(C2_WBC_eQTL_ALL$gene))
length(unique(C2_WBC_eQTL_EUR$gene))
length(unique(C2_WBC_eQTL_coloc$gene))
length(unique(C2_WBC_eQTL_coloc_eQTLgen$gene))
C2_WBC <- C2_WBC_eQTL_ALL %>% full_join(C2_WBC_eQTL_EUR, by=c("gene"="gene", "SNP"="SNP"))
C2_WBC <- C2_WBC %>% full_join(C2_WBC_eQTL_coloc, by=c("gene"="gene"))
C2_WBC <- C2_WBC %>% full_join(C2_WBC_eQTL_coloc_eQTLgen, by=c("gene"="gene"))
C2_WBC <- C2_WBC %>% mutate(outcome = "Reported infection")

WBC <- bind_rows(A2_WBC, B2_WBC, C2_WBC)
WBC <- WBC %>% mutate(OR.ALL = exp(b.x),
                      LL.ALL = exp(b.x + qnorm(0.025)*se.x),
                      UL.ALL = exp(b.x + qnorm(0.975)*se.x),
                      p.ALL = p.x,
                      OR.EUR = exp(b.y),
                      LL.EUR = exp(b.y + qnorm(0.025)*se.y),
                      UL.EUR = exp(b.y + qnorm(0.975)*se.y),
                      p.EUR = p.y,
                      NSNP.GTEx = NSNP.x,
                      PP.H0.abf.GTEx = PP.H0.abf.x, 
                      PP.H1.abf.GTEx = PP.H1.abf.x, 
                      PP.H2.abf.GTEx = PP.H2.abf.x, 
                      PP.H3.abf.GTEx = PP.H3.abf.x, 
                      PP.H4.abf.GTEx = PP.H4.abf.x,
                      NSNP.eQTLgen = NSNP.y,
                      PP.H0.abf.eQTLgen = PP.H0.abf.y,
                      PP.H1.abf.eQTLgen = PP.H1.abf.y,
                      PP.H2.abf.eQTLgen = PP.H2.abf.y,
                      PP.H3.abf.eQTLgen = PP.H3.abf.y,
                      PP.H4.abf.eQTLgen = PP.H4.abf.y
)

WBC <- WBC %>% dplyr::select(gene, outcome, SNP, OR.ALL, LL.ALL, UL.ALL, p.ALL, OR.EUR, LL.EUR, UL.EUR,
                             p.EUR, NSNP.GTEx, PP.H0.abf.GTEx, PP.H1.abf.GTEx, PP.H2.abf.GTEx, PP.H3.abf.GTEx, PP.H4.abf.GTEx,
                             NSNP.eQTLgen, PP.H0.abf.eQTLgen, PP.H1.abf.eQTLgen, PP.H2.abf.eQTLgen, PP.H3.abf.eQTLgen, PP.H4.abf.eQTLgen)

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(WBC$gene))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

WBC <- WBC %>% left_join(G_list, by=c("gene"="ensembl_gene_id"))

WBC %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable6.eQTL_WBC.xlsx")

setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

genelist <- c("ENSG00000160783", "ENSG00000160766", "ENSG00000169231", "ENSG00000185499",
              "ENSG00000168743", "ENSG00000175164", "ENSG00000068650","ENSG00000089127",
              "ENSG00000142002", "ENSG00000185303")

a_wbc_mr <- readRDS("A2_WBC_eQTL.rds")
a_wbc_mr <- a_wbc_mr %>% filter(grepl(genelist[1], exposure) | 
                                  grepl(genelist[2], exposure) | 
                                  grepl(genelist[3], exposure) | 
                                  grepl(genelist[4], exposure) |
                                  grepl(genelist[5], exposure) |
                                  grepl(genelist[6], exposure) |
                                  grepl(genelist[7], exposure) |
                                  grepl(genelist[8], exposure) |
                                  grepl(genelist[9], exposure) |
                                  grepl(genelist[10], exposure))
a_wbc_mr <- a_wbc_mr %>% drop_na(p)
a_wbc_mr <- a_wbc_mr %>% arrange(SNP)
a_wbc_mr <- a_wbc_mr %>% group_by(exposure) %>% 
  filter(b == first(b))
a_wbc_coloc <- fread("A2.coloc.wbc.eQTL.tsv")
a_wbc_mr <- a_wbc_mr %>% left_join(a_wbc_coloc, by=c("exposure"="phenotype_id"))
a_wbc_mr <- a_wbc_mr %>% mutate(outcome = "A2") %>% mutate(tissue = "WBC")


b_wbc_mr <- readRDS("B2_WBC_eQTL.rds")
b_wbc_mr <- b_wbc_mr %>% filter(grepl(genelist[1], exposure) | 
                                  grepl(genelist[2], exposure) | 
                                  grepl(genelist[3], exposure) | 
                                  grepl(genelist[4], exposure) |
                                  grepl(genelist[5], exposure) |
                                  grepl(genelist[6], exposure) |
                                  grepl(genelist[7], exposure) |
                                  grepl(genelist[8], exposure) |
                                  grepl(genelist[9], exposure) |
                                  grepl(genelist[10], exposure))
b_wbc_mr <- b_wbc_mr %>% drop_na(p)
b_wbc_mr <- b_wbc_mr %>% arrange(SNP)
b_wbc_mr <- b_wbc_mr %>% group_by(exposure) %>% 
  filter(b == first(b))
c_wbc_coloc <- fread("B2.coloc.wbc.eQTL.tsv")
b_wbc_mr <- b_wbc_mr %>% left_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
b_wbc_mr <- b_wbc_mr %>% mutate(outcome = "B2") %>% mutate(tissue = "WBC")

c_wbc_mr <- readRDS("C2_WBC_eQTL.rds")
c_wbc_mr <- c_wbc_mr %>% filter(grepl(genelist[1], exposure) | 
                             grepl(genelist[2], exposure) | 
                             grepl(genelist[3], exposure) | 
                               grepl(genelist[4], exposure) |
                               grepl(genelist[5], exposure) |
                               grepl(genelist[6], exposure) |
                               grepl(genelist[7], exposure) |
                               grepl(genelist[8], exposure) |
                               grepl(genelist[9], exposure) |
                               grepl(genelist[10], exposure))
c_wbc_mr <- c_wbc_mr %>% drop_na(p)
c_wbc_mr <- c_wbc_mr %>% arrange(SNP)
c_wbc_mr <- c_wbc_mr %>% group_by(exposure) %>% 
  filter(b == first(b))
c_wbc_coloc <- fread("C2.coloc.wbc.eQTL.tsv")
c_wbc_mr <- c_wbc_mr %>% left_join(c_wbc_coloc, by=c("exposure"="phenotype_id"))
c_wbc_mr <- c_wbc_mr %>% mutate(outcome = "C2") %>% mutate(tissue = "WBC")

a_lung_mr <- readRDS("A2_Lung_eQTL.rds")
a_lung_mr <- a_lung_mr %>% filter(grepl(genelist[1], exposure) | 
                                    grepl(genelist[2], exposure) | 
                                    grepl(genelist[3], exposure) | 
                                    grepl(genelist[4], exposure) |
                                    grepl(genelist[5], exposure) |
                                    grepl(genelist[6], exposure) |
                                    grepl(genelist[7], exposure) |
                                    grepl(genelist[8], exposure) |
                                    grepl(genelist[9], exposure) |
                                    grepl(genelist[10], exposure))

a_lung_mr <- a_lung_mr %>% drop_na(p)
a_lung_mr <- a_lung_mr %>% arrange(SNP)
a_lung_mr <- a_lung_mr %>% group_by(exposure) %>% 
  filter(b == first(b))
a_lung_coloc <- fread("A2.coloc.Lung.eQTL.tsv")
a_lung_mr <- a_lung_mr %>% left_join(c_lung_coloc, by=c("exposure"="phenotype_id"))
a_lung_mr <- a_lung_mr %>% mutate(outcome = "A2") %>% mutate(tissue = "Lung")

b_lung_mr <- readRDS("B2_Lung_eQTL.rds")
b_lung_mr <- b_lung_mr %>% filter(grepl(genelist[1], exposure) | 
                                    grepl(genelist[2], exposure) | 
                                    grepl(genelist[3], exposure) | 
                                    grepl(genelist[4], exposure) |
                                    grepl(genelist[5], exposure) |
                                    grepl(genelist[6], exposure) |
                                    grepl(genelist[7], exposure) |
                                    grepl(genelist[8], exposure) |
                                    grepl(genelist[9], exposure) |
                                    grepl(genelist[10], exposure))

b_lung_mr <- b_lung_mr %>% drop_na(p)
b_lung_mr <- b_lung_mr %>% arrange(SNP)
b_lung_mr <- b_lung_mr %>% group_by(exposure) %>% 
  filter(b == first(b))
b_lung_coloc <- fread("B2.coloc.Lung.eQTL.tsv")
b_lung_mr <- b_lung_mr %>% left_join(b_lung_coloc, by=c("exposure"="phenotype_id"))
b_lung_mr <- b_lung_mr %>% mutate(outcome = "B2") %>% mutate(tissue = "Lung")


c_lung_mr <- readRDS("C2_Lung_eQTL.rds")
c_lung_mr <- c_lung_mr %>% filter(grepl(genelist[1], exposure) | 
                                    grepl(genelist[2], exposure) | 
                                    grepl(genelist[3], exposure) | 
                                    grepl(genelist[4], exposure) |
                                    grepl(genelist[5], exposure) |
                                    grepl(genelist[6], exposure) |
                                    grepl(genelist[7], exposure) |
                                    grepl(genelist[8], exposure) |
                                    grepl(genelist[9], exposure) |
                                    grepl(genelist[10], exposure))

c_lung_mr <- c_lung_mr %>% drop_na(p)
c_lung_mr <- c_lung_mr %>% arrange(SNP)
c_lung_mr <- c_lung_mr %>% group_by(exposure) %>% 
  filter(b == first(b))
c_lung_coloc <- fread("C2.coloc.Lung.eQTL.tsv")
c_lung_mr <- c_lung_mr %>% left_join(c_lung_coloc, by=c("exposure"="phenotype_id"))
c_lung_mr <- c_lung_mr %>% mutate(outcome = "C2") %>% mutate(tissue = "Lung")
c_lung_mr %>% filter(p < 0.05 & PP.H4.abf > 0.8) %>% str()
c_mr <- bind_rows(a_wbc_mr, b_wbc_mr, c_wbc_mr, a_lung_mr, b_lung_mr, c_lung_mr)

write.xlsx(c_mr, "candidate_expression.xlsx")

c_mr <- c_mr %>% mutate(exposure = str_split(exposure, pattern="\\.", simplify = T)[,1])
c_mr <- c_mr %>% mutate(SNP = case_when(grepl("All", SNP) ~ "Inverse variance weighted",
                                        TRUE ~ SNP),
                        OR = exp(as.numeric(b)),
                        LL = exp(as.numeric(b) - qnorm(0.975)*as.numeric(se)),
                        UL = exp(as.numeric(b) + qnorm(0.975)*as.numeric(se)),
                        outcome = case_when(outcome == "A2" ~ "critical illness",
                                            outcome == "B2" ~ "hospitalization",
                                            TRUE ~ "reported SARS-CoV-2 infection"))

gene <- data.frame(unique(c_mr$exposure))
colnames(gene) <- "ensembleID"

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
colnames(gene) <- "ensembleID"


for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i],mart = ensembl)
}

c_mr <- gene %>% right_join(c_mr, by=c("ensembleID"="exposure"))

c_mr <- data.frame(c_mr)
c_mr <- c_mr %>% dplyr::select(ensembleID, gene_id, outcome, SNP, OR, LL, UL, p,
                        PP.H0.abf,    PP.H1.abf,    PP.H2.abf,   PP.H3.abf,
                        PP.H4.abf, tissue)

write.xlsx(c_mr, "/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable5.expression.xlsx")

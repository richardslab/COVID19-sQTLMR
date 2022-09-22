setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
A2_Lung_sQTL <- readRDS("A2_Lung_sQTL_ALL.rds")
A2_Lung_sQTL <- A2_Lung_sQTL %>% drop_na(p)
B2_Lung_sQTL <- readRDS("B2_Lung_sQTL_ALL.rds")
B2_Lung_sQTL <- B2_Lung_sQTL %>% drop_na(p)
C2_Lung_sQTL <- readRDS("C2_Lung_sQTL_ALL.rds")
C2_Lung_sQTL <- C2_Lung_sQTL %>% drop_na(p)

genelist <- unique(c(A2_Lung_sQTL$exposure, B2_Lung_sQTL$exposure, C2_Lung_sQTL$exposure))
length(genelist)#5999

SNPlist <- unique(c(A2_Lung_sQTL$SNP, B2_Lung_sQTL$SNP, C2_Lung_sQTL$SNP))
length(SNPlist[grepl("chr", SNPlist)])#6110

A2_Lung_sQTL <- A2_Lung_sQTL %>% arrange(SNP)
A2_Lung_sQTL_lead <- A2_Lung_sQTL[!duplicated(A2_Lung_sQTL$exposure),]
A2_Lung_sQTL_lead <- A2_Lung_sQTL_lead %>% mutate(p.adj = p.adjust(p, method = "BH"))
A2_Lung_sQTL_lead <- A2_Lung_sQTL_lead %>% dplyr::select(exposure, p.adj)

A2_Lung_sQTL <- A2_Lung_sQTL %>% dplyr::select(-p.adj) %>% left_join(A2_Lung_sQTL_lead, by="exposure")

A2_Lung_sQTL <- A2_Lung_sQTL %>% mutate(outcome = "Critical illness")

B2_Lung_sQTL <- B2_Lung_sQTL %>% arrange(SNP)
B2_Lung_sQTL_lead <- B2_Lung_sQTL[!duplicated(B2_Lung_sQTL$exposure),]
B2_Lung_sQTL_lead <- B2_Lung_sQTL_lead %>% mutate(p.adj = p.adjust(p, method = "BH"))
B2_Lung_sQTL_lead <- B2_Lung_sQTL_lead %>% dplyr::select(exposure, p.adj)

B2_Lung_sQTL <- B2_Lung_sQTL %>% dplyr::select(-p.adj) %>%  left_join(B2_Lung_sQTL_lead, by="exposure")

B2_Lung_sQTL <- B2_Lung_sQTL %>% mutate(outcome = "Hospitalization")

C2_Lung_sQTL <- C2_Lung_sQTL %>% arrange(SNP)
C2_Lung_sQTL_lead <- C2_Lung_sQTL[!duplicated(C2_Lung_sQTL$exposure),]
C2_Lung_sQTL_lead <- C2_Lung_sQTL_lead %>% mutate(p.adj = p.adjust(p, method = "BH"))
C2_Lung_sQTL_lead <- C2_Lung_sQTL_lead %>% dplyr::select(exposure, p.adj)

C2_Lung_sQTL <- C2_Lung_sQTL %>% dplyr::select(-p.adj) %>% left_join(C2_Lung_sQTL_lead, by="exposure")

C2_Lung_sQTL <- C2_Lung_sQTL %>% mutate(outcome = "Reported infection")

All_Lung_sQTL <- bind_rows(A2_Lung_sQTL, B2_Lung_sQTL, C2_Lung_sQTL)
All_Lung_sQTL <- All_Lung_sQTL %>% mutate(exposure = paste0(str_split(exposure, pattern="\\:", simplify = T)[,1],":",
                                                            str_split(exposure, pattern="\\:", simplify = T)[,2],":",
                                                            str_split(exposure, pattern="\\:", simplify = T)[,3],":",
                                                            str_split(exposure, pattern="\\:", simplify = T)[,5])) %>%
  mutate(OR = exp(b),
         LL = exp(b + qnorm(0.025)*se),
         UL = exp(b + qnorm(0.975)*se),
         ensembleID = str_split(exposure, pattern="\\:", simplify = T)[,4],
         ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(All_Lung_sQTL$ensembleID))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

All_Lung_sQTL <- All_Lung_sQTL %>% left_join(G_list, by=c("ensembleID"="ensembl_gene_id"))

saveRDS(All_Lung_sQTL, "All_Lung_sQTL_ALL.rds")

All_Lung_sQTL <- readRDS("All_Lung_sQTL_ALL.rds")

length(unique(All_Lung_sQTL$exposure))
length(unique(All_Lung_sQTL$SNP)) - 2
length(unique(All_Lung_sQTL$ensembleID))

All_Lung_sQTL %>% filter(p < 0.05/9768) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_Lung_sQTL %>% filter(p < 0.05/9768) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_Lung_sQTL %>% filter(p < 0.05/9768) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()
18+22+20

All_Lung_sQTL %>% filter(p < 0.05/9768) %>% dplyr::select(ensembleID) %>% unique() %>% dim()



final <- All_Lung_sQTL %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, SNP, OR, LL, UL, p, p.adj)
final %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable1.Lung.MR.xlsx")


A2_WBC_sQTL <- readRDS("A2_WBC_sQTL_ALL.rds")
A2_WBC_sQTL <- A2_WBC_sQTL %>% drop_na(p)
B2_WBC_sQTL <- readRDS("B2_WBC_sQTL_ALL.rds")
B2_WBC_sQTL <- B2_WBC_sQTL %>% drop_na(p)
C2_WBC_sQTL <- readRDS("C2_WBC_sQTL_ALL.rds")
C2_WBC_sQTL <- C2_WBC_sQTL %>% drop_na(p)

genelist <- unique(c(A2_WBC_sQTL$exposure, B2_WBC_sQTL$exposure, C2_WBC_sQTL$exposure))
length(genelist)#3769

SNPlist <- unique(c(A2_WBC_sQTL$SNP, B2_WBC_sQTL$SNP, C2_WBC_sQTL$SNP))
length(SNPlist[grepl("chr", SNPlist)])#3883


A2_WBC_sQTL <- A2_WBC_sQTL %>% arrange(SNP)
A2_WBC_sQTL_lead <- A2_WBC_sQTL[!duplicated(A2_WBC_sQTL$exposure),]
A2_WBC_sQTL_lead <- A2_WBC_sQTL_lead %>% mutate(p.adj = p.adjust(p, method = "BH"))
A2_WBC_sQTL_lead <- A2_WBC_sQTL_lead %>% dplyr::select(exposure, p.adj)

A2_WBC_sQTL <- A2_WBC_sQTL %>% left_join(A2_WBC_sQTL_lead, by="exposure")

A2_WBC_sQTL <- A2_WBC_sQTL %>% mutate(outcome = "Critical illness")

B2_WBC_sQTL <- B2_WBC_sQTL %>% arrange(SNP)
B2_WBC_sQTL_lead <- B2_WBC_sQTL[!duplicated(B2_WBC_sQTL$exposure),]
B2_WBC_sQTL_lead <- B2_WBC_sQTL_lead %>% mutate(p.adj = p.adjust(p, method = "BH"))
B2_WBC_sQTL_lead <- B2_WBC_sQTL_lead %>% dplyr::select(exposure, p.adj)

B2_WBC_sQTL <- B2_WBC_sQTL  %>%  dplyr::select(-p.adj) %>% left_join(B2_WBC_sQTL_lead, by="exposure")

B2_WBC_sQTL <- B2_WBC_sQTL %>% mutate(outcome = "Hospitalization")

C2_WBC_sQTL <- C2_WBC_sQTL %>% arrange(SNP)
C2_WBC_sQTL_lead <- C2_WBC_sQTL[!duplicated(C2_WBC_sQTL$exposure),]
C2_WBC_sQTL_lead <- C2_WBC_sQTL_lead %>% mutate(p.adj = p.adjust(p, method = "BH"))
C2_WBC_sQTL_lead <- C2_WBC_sQTL_lead %>% dplyr::select(exposure, p.adj)

C2_WBC_sQTL <- C2_WBC_sQTL %>%  dplyr::select(-p.adj) %>%  left_join(C2_WBC_sQTL_lead, by="exposure")

C2_WBC_sQTL <- C2_WBC_sQTL %>% mutate(outcome = "Reported infection")

All_WBC_sQTL <- bind_rows(A2_WBC_sQTL, B2_WBC_sQTL, C2_WBC_sQTL)
All_WBC_sQTL <- All_WBC_sQTL %>% mutate(exposure = paste0(str_split(exposure, pattern="\\:", simplify = T)[,1],":",
                                                            str_split(exposure, pattern="\\:", simplify = T)[,2],":",
                                                            str_split(exposure, pattern="\\:", simplify = T)[,3],":",
                                                            str_split(exposure, pattern="\\:", simplify = T)[,5])) %>%
  mutate(OR = exp(b),
         LL = exp(b + qnorm(0.025)*se),
         UL = exp(b + qnorm(0.975)*se),
         ensembleID = str_split(exposure, pattern="\\:", simplify = T)[,4],
         ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(All_WBC_sQTL$ensembleID))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

All_WBC_sQTL <- All_WBC_sQTL %>% left_join(G_list, by=c("ensembleID"="ensembl_gene_id"))

saveRDS(All_WBC_sQTL, "All_WBC_sQTL_ALL.rds")

All_WBC_sQTL %>% filter(p < 0.05/9768) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_WBC_sQTL %>% filter(p < 0.05/9768) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_WBC_sQTL %>% filter(p < 0.05/9768) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()
12+13+9

All_WBC_sQTL %>% filter(p < 0.05/9768)  %>% dplyr::select(ensembleID) %>% unique() %>% dim()


length(unique(All_WBC_sQTL$exposure))
length(unique(All_WBC_sQTL$SNP)) - 2
length(unique(All_WBC_sQTL$ensembleID))

final <- All_WBC_sQTL %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, SNP, OR, LL, UL, p, p.adj)
final %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable2.WBC.MR.xlsx")

All_Lung_sQTL <- readRDS("All_Lung_sQTL_ALL.rds")
All_WBC_sQTL <- readRDS("All_WBC_sQTL_ALL.rds")

ALL <- bind_rows(All_WBC_sQTL, All_Lung_sQTL)

length(unique(ALL$exposure))
length(unique(ALL$SNP)) - 2
length(unique(ALL$ensembleID))



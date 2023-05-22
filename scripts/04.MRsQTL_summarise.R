setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
A2_Lung_sQTL <- readRDS("A2_Lung_sQTL_EUR.rds")
A2_Lung_sQTL <- A2_Lung_sQTL %>% drop_na(p)
length(unique(A2_Lung_sQTL$exposure))#5546
B2_Lung_sQTL <- readRDS("B2_Lung_sQTL_EUR.rds")
B2_Lung_sQTL <- B2_Lung_sQTL %>% drop_na(p)
length(unique(B2_Lung_sQTL$exposure))#5570
C2_Lung_sQTL <- readRDS("C2_Lung_sQTL_EUR.rds")
C2_Lung_sQTL <- C2_Lung_sQTL %>% drop_na(p)
length(unique(C2_Lung_sQTL$exposure))#5686

genelist <- unique(c(A2_Lung_sQTL$exposure, B2_Lung_sQTL$exposure, C2_Lung_sQTL$exposure))
length(genelist)#5724

SNPlist <- unique(c(A2_Lung_sQTL$SNP, B2_Lung_sQTL$SNP, C2_Lung_sQTL$SNP))
length(SNPlist[grepl("chr", SNPlist)])#5807

A2_Lung_sQTL <- A2_Lung_sQTL %>% arrange(SNP)
A2_Lung_sQTL <- A2_Lung_sQTL %>% mutate(outcome = "Critical illness")

A2_Lung_sQTL <- A2_Lung_sQTL %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=","))
A2_Lung_sQTL <- A2_Lung_sQTL %>% 
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
A2_Lung_sQTL <- A2_Lung_sQTL %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))
A2_Lung_sQTL <- A2_Lung_sQTL %>% group_by(exposure) %>% 
  filter(Method == first(Method))


B2_Lung_sQTL <- B2_Lung_sQTL %>% arrange(SNP)
B2_Lung_sQTL <- B2_Lung_sQTL %>% mutate(outcome = "Hospitalization")

B2_Lung_sQTL <- B2_Lung_sQTL %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=","))
B2_Lung_sQTL <- B2_Lung_sQTL %>% 
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
B2_Lung_sQTL <- B2_Lung_sQTL %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))
B2_Lung_sQTL <- B2_Lung_sQTL %>% group_by(exposure) %>% 
  filter(Method == first(Method))

C2_Lung_sQTL <- C2_Lung_sQTL %>% arrange(SNP)
C2_Lung_sQTL <- C2_Lung_sQTL %>% mutate(outcome = "Reported infection")

C2_Lung_sQTL <- C2_Lung_sQTL %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=","))
C2_Lung_sQTL <- C2_Lung_sQTL %>% 
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
C2_Lung_sQTL <- C2_Lung_sQTL %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))
C2_Lung_sQTL <- C2_Lung_sQTL %>% group_by(exposure) %>% 
  filter(Method == first(Method))


All_Lung_sQTL <- bind_rows(A2_Lung_sQTL, B2_Lung_sQTL, C2_Lung_sQTL)
All_Lung_sQTL <- All_Lung_sQTL %>% mutate(OR = exp(b),
         LL = exp(b + qnorm(0.025)*se),
         UL = exp(b + qnorm(0.975)*se),
         ensembleID = str_split(exposure, pattern="\\:", simplify = T)[,5],
         ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])
#All_Lung_sQTL <- All_Lung_sQTL %>% mutate(ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(All_Lung_sQTL$ensembleID))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

All_Lung_sQTL <- All_Lung_sQTL %>% left_join(G_list, by=c("ensembleID"="ensembl_gene_id"))

saveRDS(All_Lung_sQTL, "All_Lung_sQTL_EUR.rds")

All_Lung_sQTL <- readRDS("All_Lung_sQTL_EUR.rds")

length(unique(All_Lung_sQTL$exposure))
length(unique(All_Lung_sQTL$SNP)) - 2
length(unique(All_Lung_sQTL$ensembleID))

All_Lung_sQTL %>% filter(p < 0.05/27230) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_Lung_sQTL %>% filter(p < 0.05/27230) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_Lung_sQTL %>% filter(p < 0.05/27230) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()
 

All_Lung_sQTL %>% filter(p < 0.05/27230) %>% dplyr::select(ensembleID) %>% unique() %>% dim()


final <- All_Lung_sQTL %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, Method, SNPlist, OR, LL, UL, p)
final %>% openxlsx::write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable1.Lung.MR.xlsx")


A2_WBC_sQTL <- readRDS("A2_WBC_sQTL_EUR.rds")
A2_WBC_sQTL <- A2_WBC_sQTL %>% drop_na(p)
length(unique(A2_WBC_sQTL$exposure))#3453
B2_WBC_sQTL <- readRDS("B2_WBC_sQTL_EUR.rds")
B2_WBC_sQTL <- B2_WBC_sQTL %>% drop_na(p)
length(unique(B2_WBC_sQTL$exposure))#3451
C2_WBC_sQTL <- readRDS("C2_WBC_sQTL_EUR.rds")
C2_WBC_sQTL <- C2_WBC_sQTL %>% drop_na(p)
length(unique(C2_WBC_sQTL$exposure))#3524

5546+5570+5686+3453+3451+3524#27230


genelist <- unique(c(A2_Lung_sQTL$exposure, B2_Lung_sQTL$exposure, C2_Lung_sQTL$exposure, A2_WBC_sQTL$exposure, B2_WBC_sQTL$exposure, C2_WBC_sQTL$exposure))
length(genelist)#9292

SNPlist <- unique(c(A2_Lung_sQTL$SNP, B2_Lung_sQTL$SNP, C2_Lung_sQTL$SNP, A2_WBC_sQTL$SNP, B2_WBC_sQTL$SNP, C2_WBC_sQTL$SNP))
length(SNPlist[grepl("chr", SNPlist)])#3658

SNPlist <- unique(c(A2_Lung_sQTL$SNP, B2_Lung_sQTL$SNP, C2_Lung_sQTL$SNP, A2_WBC_sQTL$SNP, B2_WBC_sQTL$SNP, C2_WBC_sQTL$SNP))
length(SNPlist[grepl("chr", SNPlist)])#3658


A2_WBC_sQTL <- A2_WBC_sQTL %>% arrange(SNP)
A2_WBC_sQTL <- A2_WBC_sQTL %>% mutate(outcome = "Critical illness")

A2_WBC_sQTL <- A2_WBC_sQTL %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=","))
A2_WBC_sQTL <- A2_WBC_sQTL %>% 
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
A2_WBC_sQTL <- A2_WBC_sQTL %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))
A2_WBC_sQTL <- A2_WBC_sQTL %>% group_by(exposure) %>% 
  filter(Method == first(Method))


B2_WBC_sQTL <- B2_WBC_sQTL %>% arrange(SNP)
B2_WBC_sQTL <- B2_WBC_sQTL %>% mutate(outcome = "Hospitalization")

B2_WBC_sQTL <- B2_WBC_sQTL %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=","))
B2_WBC_sQTL <- B2_WBC_sQTL %>% 
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
B2_WBC_sQTL <- B2_WBC_sQTL %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))
B2_WBC_sQTL <- B2_WBC_sQTL %>% group_by(exposure) %>% 
  filter(Method == first(Method))

C2_WBC_sQTL <- C2_WBC_sQTL %>% arrange(SNP)
C2_WBC_sQTL <- C2_WBC_sQTL %>% mutate(outcome = "Reported infection")

C2_WBC_sQTL <- C2_WBC_sQTL %>% group_by(exposure) %>% 
  mutate(SNPlist = paste0(SNP, collapse=","))
C2_WBC_sQTL <- C2_WBC_sQTL %>% 
  mutate(SNPlist = gsub("All - Inverse variance weighted,","", SNPlist),
         SNPlist = gsub("All - MR Egger,","", SNPlist)) %>% ungroup() 
C2_WBC_sQTL <- C2_WBC_sQTL %>% mutate(Method = ifelse(grepl("chr", SNP), "Wald ratio", "Inverse variance weighted"))
C2_WBC_sQTL <- C2_WBC_sQTL %>% group_by(exposure) %>% 
  filter(Method == first(Method))


All_WBC_sQTL <- bind_rows(A2_WBC_sQTL, B2_WBC_sQTL, C2_WBC_sQTL)
All_WBC_sQTL <- All_WBC_sQTL  %>%
  mutate(OR = exp(b),
         LL = exp(b + qnorm(0.025)*se),
         UL = exp(b + qnorm(0.975)*se),
         ensembleID = str_split(exposure, pattern="\\:", simplify = T)[,4],
         ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])
All_WBC_sQTL <- All_WBC_sQTL  %>%
  mutate(ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])

tmp <- readRDS("Whole_Blood_cluster_intron_gene_map.rds")
tmp <- tmp %>% mutate(exposure = paste0(pos,":",gene),
                      exposure1 = paste0(pos,":",cluster,":",gene)) %>% dplyr::select(exposure, exposure1)
All_WBC_sQTL <- All_WBC_sQTL %>% inner_join(tmp, by="exposure")
All_WBC_sQTL <- All_WBC_sQTL %>% mutate(exposure = exposure1)


library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(All_WBC_sQTL$ensembleID))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

All_WBC_sQTL <- All_WBC_sQTL %>% left_join(G_list, by=c("ensembleID"="ensembl_gene_id"))

saveRDS(All_WBC_sQTL, "All_WBC_sQTL_EUR.rds")

All_WBC_sQTL %>% filter(p < 0.05/27230) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_WBC_sQTL %>% filter(p < 0.05/27230) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
All_WBC_sQTL %>% filter(p < 0.05/27230) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()
7+12+7

All_WBC_sQTL %>% filter(p < 0.05/27230)  %>% dplyr::select(ensembleID) %>% unique() %>% dim()


length(unique(All_WBC_sQTL$exposure))
length(unique(All_WBC_sQTL$SNP)) - 2
length(unique(c(All_WBC_sQTL$ensembleID, All_Lung_sQTL$ensembleID)))#5097

final <- All_WBC_sQTL %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, Method, SNPlist, OR, LL, UL, p, p.adj)
final %>% openxlsx::write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable2.WBC.MR.xlsx")




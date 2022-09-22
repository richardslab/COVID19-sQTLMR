setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

lung <- read.xlsx("Lung.coloc.xlsx")
lunglist <- lung %>% filter(PP.H4.abf > PP.H0.abf + PP.H1.abf + PP.H2.abf + PP.H3.abf)
lunglist <- unique(lunglist$junction)

sQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Lung.v8.independent_sqtls.txt.gz")
sQTL_Lung <- sQTL_Lung %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                  POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                  EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                  NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
sQTL_Lung <- sQTL_Lung %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                  N = round(ma_count/maf*(1/2)))

#sQTL_Lung %>% filter(grepl("chr13:112875941:112880546:clu_3196:ENSG00000068650.18", phenotype_id))
mrresult <- read.xlsx("All_Lung_sQTL.xlsx")
mrresult_a <- mrresult %>% filter(outcome == "A2")#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_a <- mrresult_a %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_a <- unique(mrresult_a) %>% drop_na(p)
a2 <- fread("A2.coloc.lung.tsv")
a2 <- a2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
a2 <- unique(a2)
mrresult_a <- mrresult_a %>% left_join(a2, by="junction")
mrresult_a <- mrresult_a %>% filter(junction %in% lunglist)

mrresult_b <- mrresult %>% filter(outcome == "B2") #  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_b <- mrresult_b %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_b <- unique(mrresult_b) %>% drop_na(p)

b2 <- fread("B2.coloc.lung.tsv")
b2 <- b2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
mrresult_b <- mrresult_b %>% left_join(b2, by="junction")
mrresult_b <- mrresult_b %>% filter(junction %in% lunglist)

mrresult_c <- mrresult %>% filter(outcome == "C2")# %>% filter(p.adj < 0.05)#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_c <- mrresult_c %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_c <- unique(mrresult_c) %>% drop_na(p)
c2 <- fread("C2.coloc.lung.tsv")
c2 <- c2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
mrresult_c <- mrresult_c %>% left_join(c2, by="junction")
mrresult_c <- mrresult_c %>% filter(junction %in% lunglist)

mrresult <- bind_rows(mrresult_a, mrresult_b, mrresult_c)

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene <- data.frame(unique(mrresult$gene))
colnames(gene) <- "ensembleID"

for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i], mart = ensembl)
}

mrresult <- gene %>% right_join(mrresult, by=c("ensembleID"="gene"))

mrresult <- mrresult %>% mutate(splicing = paste0(gene_id," (", str_split(junction, pattern=":", simplify = T)[,1], ":",
                                                  str_split(junction, pattern=":", simplify = T)[,2], "-",
                                                  str_split(junction, pattern=":", simplify = T)[,3],")"))

mrresult <- mrresult %>% dplyr::select(splicing, ensembleID, outcome, SNP, OR, LL, UL, p, PP.H0.abf, PP.H1.abf,
                                       PP.H2.abf, PP.H3.abf, PP.H4.abf)

mrresult_lung <- mrresult %>% mutate(tissue = "Lung")

wbc <- read.xlsx("WBC.coloc.xlsx")
wbclist <- wbc %>% filter(PP.H4.abf > PP.H0.abf + PP.H1.abf + PP.H2.abf + PP.H3.abf)
wbclist <- unique(wbclist$junction)

mrresult <- read.xlsx("All_WBC_sQTL.xlsx")
mrresult_a <- mrresult %>% filter(outcome == "A2")#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_a <- mrresult_a %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_a <- unique(mrresult_a) %>% drop_na(p)
a2 <- fread("A2.coloc.lung.tsv")
a2 <- a2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
a2 <- unique(a2)
mrresult_a <- mrresult_a %>% left_join(a2, by="junction")
mrresult_a <- mrresult_a %>% filter(junction %in% wbclist)

mrresult_b <- mrresult %>% filter(outcome == "B2") #  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_b <- mrresult_b %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_b <- unique(mrresult_b) %>% drop_na(p)

b2 <- fread("B2.coloc.lung.tsv")
b2 <- b2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
mrresult_b <- mrresult_b %>% left_join(b2, by="junction")
mrresult_b <- mrresult_b %>% filter(junction %in% wbclist)

mrresult_c <- mrresult %>% filter(outcome == "C2")# %>% filter(p.adj < 0.05)#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_c <- mrresult_c %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_c <- unique(mrresult_c) %>% drop_na(p)
c2 <- fread("C2.coloc.lung.tsv")
c2 <- c2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
mrresult_c <- mrresult_c %>% left_join(c2, by="junction")
mrresult_c <- mrresult_c %>% filter(junction %in% wbclist)

mrresult <- bind_rows(mrresult_a, mrresult_b, mrresult_c)

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene <- data.frame(unique(mrresult$gene))
colnames(gene) <- "ensembleID"

for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i], mart = ensembl)
}

mrresult <- gene %>% right_join(mrresult, by=c("ensembleID"="gene"))

mrresult <- mrresult %>% mutate(splicing = paste0(gene_id," (", str_split(junction, pattern=":", simplify = T)[,1], ":",
                                                  str_split(junction, pattern=":", simplify = T)[,2], "-",
                                                  str_split(junction, pattern=":", simplify = T)[,3],")"))

mrresult <- mrresult %>% dplyr::select(splicing, ensembleID, outcome, SNP, OR, LL, UL, p, PP.H0.abf, PP.H1.abf,
                                       PP.H2.abf, PP.H3.abf, PP.H4.abf)

mrresult_wbc <- mrresult %>% mutate(tissue = "WBC")


mrresult <- bind_rows(mrresult_lung, mrresult_wbc)

write.xlsx(mrresult, "significant.MR.coloc.xlsx")
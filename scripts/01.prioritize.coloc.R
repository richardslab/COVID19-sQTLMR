setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/04.Barreiro/coloc_results/")
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(coloc)

results <- fread("high.coloc.summary.tsv.gz")
results <- results %>% arrange(gene_id)
write.xlsx(results, file="high.coloc.xlsx")

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(results$gene_id))
colnames(gene) <- "ensembleID"
G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'hgnc_symbol',values = gene$ensembleID,mart = ensembl)

GTEx_A2 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GTEx/eQTL/A2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz")
GTEx_A2 <- GTEx_A2 %>% mutate(geneid = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
GTEx_A2 <- GTEx_A2 %>% filter(geneid %in% G_list$ensembl_gene_id) %>% mutate(pheno = "A2_EUR")
GTEx_B2 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GTEx/eQTL/B2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz")
GTEx_B2 <- GTEx_B2 %>% mutate(geneid = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
GTEx_B2 <- GTEx_B2 %>% filter(geneid %in% G_list$ensembl_gene_id) %>% mutate(pheno = "B2_EUR")
GTEx_C2 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GTEx/eQTL/C2_eQTL_GTEx_Whole_Blood_coloc.tsv.gz")
GTEx_C2 <- GTEx_C2 %>% mutate(geneid = str_split(phenotype_id, pattern="\\.", simplify = T)[,1])
GTEx_C2 <- GTEx_C2 %>% filter(geneid %in% G_list$ensembl_gene_id) %>% mutate(pheno = "C2_EUR")

GTEx <- bind_rows(GTEx_A2,GTEx_B2, GTEx_C2)
GTEx <- GTEx %>% left_join(G_list, by=c("geneid"="ensembl_gene_id"))

write.xlsx(GTEx, file="GTEx.coloc.xlsx")

eQTLgen_A2 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/eQTLgen/A2_eQTL_eQTLgen.tsv")
eQTLgen_A2 <- eQTLgen_A2 %>% filter(phenotype_id %in% G_list$ensembl_gene_id) %>% mutate(pheno = "A2_EUR")
eQTLgen_B2 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/eQTLgen/B2_eQTL_eQTLgen.tsv")
eQTLgen_B2 <- eQTLgen_B2 %>% filter(phenotype_id %in% G_list$ensembl_gene_id) %>% mutate(pheno = "B2_EUR")
eQTLgen_C2 <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/eQTLgen/C2_eQTL_eQTLgen.tsv")
eQTLgen_C2 <- eQTLgen_C2 %>% filter(phenotype_id %in% G_list$ensembl_gene_id) %>% mutate(pheno = "C2_EUR")

eQTLgen <- bind_rows(eQTLgen_A2,eQTLgen_B2, eQTLgen_C2)
eQTLgen <- eQTLgen %>% left_join(G_list, by=c("phenotype_id"="ensembl_gene_id"))

write.xlsx(eQTLgen, file="eQTLgen.coloc.xlsx")




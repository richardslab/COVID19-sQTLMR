setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

#sQTL
sQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Lung.v8.independent_sqtls.txt.gz")

tmp <- sQTL_Lung %>% select(phenotype_id)
tmp <- unique(tmp)
tmp <- tmp %>% mutate(pos = paste0(str_split(phenotype_id, pattern=":", simplify = T)[,1],":",
                                   str_split(phenotype_id, pattern=":", simplify = T)[,2],":",
                                   str_split(phenotype_id, pattern=":", simplify = T)[,3]),
                      cluster = str_split(phenotype_id, pattern=":", simplify = T)[,4],
                      gene = str_split(phenotype_id, pattern=":", simplify = T)[,5])

tmp <- tmp %>% select(pos, cluster, gene)
saveRDS(tmp, file="Lung_cluster_intron_gene_map.rds")                  
                      
# sQTL_Lung_EUR_chr1 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr1.tsv.gz")
# sQTL_Lung_EUR_chr1 <- sQTL_Lung_EUR_chr1 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr1 <- sQTL_Lung_EUR_chr1 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))

sQTL_Lung <- sQTL_Lung %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
                                                        str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
                                                        str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
                                                        str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))

# sQTL_Lung_EUR_chr1 <- sQTL_Lung_EUR_chr1 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr1 <- sQTL_Lung_EUR_chr1 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr1 <- sQTL_Lung_EUR_chr1 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr2 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr2.tsv.gz")
# sQTL_Lung_EUR_chr2 <- sQTL_Lung_EUR_chr2 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr2 <- sQTL_Lung_EUR_chr2 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr2 <- sQTL_Lung_EUR_chr2 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr2 <- sQTL_Lung_EUR_chr2 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr2 <- sQTL_Lung_EUR_chr2 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr3 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr3.tsv.gz")
# sQTL_Lung_EUR_chr3 <- sQTL_Lung_EUR_chr3 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr3 <- sQTL_Lung_EUR_chr3 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr3 <- sQTL_Lung_EUR_chr3 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr3 <- sQTL_Lung_EUR_chr3 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr3 <- sQTL_Lung_EUR_chr3 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr4 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr4.tsv.gz")
# sQTL_Lung_EUR_chr4 <- sQTL_Lung_EUR_chr4 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr4 <- sQTL_Lung_EUR_chr4 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr4 <- sQTL_Lung_EUR_chr4 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr4 <- sQTL_Lung_EUR_chr4 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr4 <- sQTL_Lung_EUR_chr4 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr5 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr5.tsv.gz")
# sQTL_Lung_EUR_chr5 <- sQTL_Lung_EUR_chr5 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr5 <- sQTL_Lung_EUR_chr5 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr5 <- sQTL_Lung_EUR_chr5 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr5 <- sQTL_Lung_EUR_chr5 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr5 <- sQTL_Lung_EUR_chr5 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr6 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr6.tsv.gz")
# sQTL_Lung_EUR_chr6 <- sQTL_Lung_EUR_chr6 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr6 <- sQTL_Lung_EUR_chr6 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr6 <- sQTL_Lung_EUR_chr6 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr6 <- sQTL_Lung_EUR_chr6 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr6 <- sQTL_Lung_EUR_chr6 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr7 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr7.tsv.gz")
# sQTL_Lung_EUR_chr7 <- sQTL_Lung_EUR_chr7 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr7 <- sQTL_Lung_EUR_chr7 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr7 <- sQTL_Lung_EUR_chr7 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr7 <- sQTL_Lung_EUR_chr7 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr7 <- sQTL_Lung_EUR_chr7 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr8 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr8.tsv.gz")
# sQTL_Lung_EUR_chr8 <- sQTL_Lung_EUR_chr8 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr8 <- sQTL_Lung_EUR_chr8 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr8 <- sQTL_Lung_EUR_chr8 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr8 <- sQTL_Lung_EUR_chr8 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr8 <- sQTL_Lung_EUR_chr8 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr9 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr9.tsv.gz")
# sQTL_Lung_EUR_chr9 <- sQTL_Lung_EUR_chr9 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr9 <- sQTL_Lung_EUR_chr9 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr9 <- sQTL_Lung_EUR_chr9 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr9 <- sQTL_Lung_EUR_chr9 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr9 <- sQTL_Lung_EUR_chr9 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr10 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr10.tsv.gz")
# sQTL_Lung_EUR_chr10 <- sQTL_Lung_EUR_chr10 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr10 <- sQTL_Lung_EUR_chr10 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr10 <- sQTL_Lung_EUR_chr10 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr10 <- sQTL_Lung_EUR_chr10 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr10 <- sQTL_Lung_EUR_chr10 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr11 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr11.tsv.gz")
# sQTL_Lung_EUR_chr11 <- sQTL_Lung_EUR_chr11 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr11 <- sQTL_Lung_EUR_chr11 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr11 <- sQTL_Lung_EUR_chr11 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr11 <- sQTL_Lung_EUR_chr11 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr11 <- sQTL_Lung_EUR_chr11 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr12 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr12.tsv.gz")
# sQTL_Lung_EUR_chr12 <- sQTL_Lung_EUR_chr12 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr12 <- sQTL_Lung_EUR_chr12 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr12 <- sQTL_Lung_EUR_chr12 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr12 <- sQTL_Lung_EUR_chr12 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr12 <- sQTL_Lung_EUR_chr12 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr13 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr13.tsv.gz")
# sQTL_Lung_EUR_chr13 <- sQTL_Lung_EUR_chr13 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr13 <- sQTL_Lung_EUR_chr13 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr13 <- sQTL_Lung_EUR_chr13 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr13 <- sQTL_Lung_EUR_chr13 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr13 <- sQTL_Lung_EUR_chr13 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr14 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr14.tsv.gz")
# sQTL_Lung_EUR_chr14 <- sQTL_Lung_EUR_chr14 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr14 <- sQTL_Lung_EUR_chr14 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr14 <- sQTL_Lung_EUR_chr14 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr14 <- sQTL_Lung_EUR_chr14 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr14 <- sQTL_Lung_EUR_chr14 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr15 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr15.tsv.gz")
# sQTL_Lung_EUR_chr15 <- sQTL_Lung_EUR_chr15 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr15 <- sQTL_Lung_EUR_chr15 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr15 <- sQTL_Lung_EUR_chr15 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr15 <- sQTL_Lung_EUR_chr15 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr15 <- sQTL_Lung_EUR_chr15 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr16 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr16.tsv.gz")
# sQTL_Lung_EUR_chr16 <- sQTL_Lung_EUR_chr16 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr16 <- sQTL_Lung_EUR_chr16 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr16 <- sQTL_Lung_EUR_chr16 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr16 <- sQTL_Lung_EUR_chr16 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr16 <- sQTL_Lung_EUR_chr16 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr17 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr17.tsv.gz")
# sQTL_Lung_EUR_chr17 <- sQTL_Lung_EUR_chr17 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr17 <- sQTL_Lung_EUR_chr17 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr17 <- sQTL_Lung_EUR_chr17 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr17 <- sQTL_Lung_EUR_chr17 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr17 <- sQTL_Lung_EUR_chr17 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr18 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr18.tsv.gz")
# sQTL_Lung_EUR_chr18 <- sQTL_Lung_EUR_chr18 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr18 <- sQTL_Lung_EUR_chr18 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr18 <- sQTL_Lung_EUR_chr18 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr18 <- sQTL_Lung_EUR_chr18 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr18 <- sQTL_Lung_EUR_chr18 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr19 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr19.tsv.gz")
# sQTL_Lung_EUR_chr19 <- sQTL_Lung_EUR_chr19 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr19 <- sQTL_Lung_EUR_chr19 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr19 <- sQTL_Lung_EUR_chr19 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr19 <- sQTL_Lung_EUR_chr19 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr19 <- sQTL_Lung_EUR_chr19 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr20 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr20.tsv.gz")
# sQTL_Lung_EUR_chr20 <- sQTL_Lung_EUR_chr20 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr20 <- sQTL_Lung_EUR_chr20 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr20 <- sQTL_Lung_EUR_chr20 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr20 <- sQTL_Lung_EUR_chr20 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr20 <- sQTL_Lung_EUR_chr20 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr21 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr21.tsv.gz")
# sQTL_Lung_EUR_chr21 <- sQTL_Lung_EUR_chr21 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr21 <- sQTL_Lung_EUR_chr21 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr21 <- sQTL_Lung_EUR_chr21 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr21 <- sQTL_Lung_EUR_chr21 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr21 <- sQTL_Lung_EUR_chr21 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR_chr22 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr22.tsv.gz")
# sQTL_Lung_EUR_chr22 <- sQTL_Lung_EUR_chr22 %>% filter(variant_id %in% sQTL_Lung$variant_id)
# sQTL_Lung_EUR_chr22 <- sQTL_Lung_EUR_chr22 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_Lung_EUR_chr22 <- sQTL_Lung_EUR_chr22 %>% filter(phenotype_id %in% sQTL_Lung$phenotype_id)
#
# sQTL_Lung_EUR_chr22 <- sQTL_Lung_EUR_chr22 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_Lung_EUR_chr22 <- sQTL_Lung_EUR_chr22 %>% inner_join(sQTL_Lung, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_Lung_EUR <- bind_rows(sQTL_Lung_EUR_chr1, sQTL_Lung_EUR_chr2, sQTL_Lung_EUR_chr3, sQTL_Lung_EUR_chr4, sQTL_Lung_EUR_chr5,
#                            sQTL_Lung_EUR_chr6, sQTL_Lung_EUR_chr7, sQTL_Lung_EUR_chr8, sQTL_Lung_EUR_chr9, sQTL_Lung_EUR_chr10,
#                            sQTL_Lung_EUR_chr11, sQTL_Lung_EUR_chr12, sQTL_Lung_EUR_chr13, sQTL_Lung_EUR_chr14, sQTL_Lung_EUR_chr15,
#                            sQTL_Lung_EUR_chr16, sQTL_Lung_EUR_chr17, sQTL_Lung_EUR_chr18, sQTL_Lung_EUR_chr19, sQTL_Lung_EUR_chr20,
#                            sQTL_Lung_EUR_chr21, sQTL_Lung_EUR_chr22)
#
# sQTL_Lung_EUR <- sQTL_Lung_EUR %>% mutate(N = round(ma_count.x/maf.x*(1/2)))
# sQTL_Lung_EUR <- sQTL_Lung_EUR %>% mutate(SNP = paste0(CHR,":",POS))
# exp_sQTL_Lung_EUR <- format_data(sQTL_Lung_EUR, type="exposure",
#                                  phenotype_col = "phenotype_id",
#                                  snp_col = "SNP",
#                                  beta_col = "slope.x",
#                                  se_col = "slope_se.x",
#                                  effect_allele_col = "EA",
#                                  other_allele_col = "NEA",
#                                  eaf = "maf.x",
#                                  pval_col = "pval_nominal.x",
#                                  samplesize_col = "N",
#                                  chr_col = "CHR",
#                                  pos_col = "POS",
# )
exp_sQTL_Lung_EUR <- exp_sQTL_Lung_EUR %>% filter(!(chr.exposure == "chr6" & pos.exposure > 28510120 & pos.exposure < 33480577))

# saveRDS(exp_sQTL_Lung_EUR, file="exposure_sQTL_Lung_EUR.rds")
# exp_sQTL_Lung_EUR <- readRDS(file="exposure_sQTL_Lung_EUR.rds")
# 
# tmp <- readRDS(file="Lung_cluster_intron_gene_map.rds")  
# tmp <- tmp %>% mutate(exposure = paste0(pos,":",gene),
#                       exposure1 = paste0(pos,":",cluster,":",gene)) %>% select(exposure, exposure1)
# exp_sQTL_Lung_EUR <- exp_sQTL_Lung_EUR %>% inner_join(tmp, by="exposure")
# exp_sQTL_Lung_EUR <- exp_sQTL_Lung_EUR %>% mutate(exposure = exposure1)

saveRDS(exp_sQTL_Lung_EUR, file="exposure_sQTL_Lung_EUR.rds")

sQTL_WBC <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Whole_Blood.v8.independent_sqtls.txt.gz")

tmp <- sQTL_WBC %>% select(phenotype_id)
tmp <- unique(tmp)
tmp <- tmp %>% mutate(pos = paste0(str_split(phenotype_id, pattern=":", simplify = T)[,1],":",
                                   str_split(phenotype_id, pattern=":", simplify = T)[,2],":",
                                   str_split(phenotype_id, pattern=":", simplify = T)[,3]),
                      cluster = str_split(phenotype_id, pattern=":", simplify = T)[,4],
                      gene = str_split(phenotype_id, pattern=":", simplify = T)[,5])
tmp <- tmp %>% select(pos, cluster, gene)
saveRDS(tmp, file="Whole_Blood_cluster_intron_gene_map.rds")    

# sQTL_WBC_EUR_chr1 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr1.tsv.gz")
# sQTL_WBC_EUR_chr1 <- sQTL_WBC_EUR_chr1 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr1 <- sQTL_WBC_EUR_chr1 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))

sQTL_WBC <- sQTL_WBC %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
                                                        str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
                                                        str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
                                                        str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))

# sQTL_WBC_EUR_chr1 <- sQTL_WBC_EUR_chr1 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr1 <- sQTL_WBC_EUR_chr1 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr1 <- sQTL_WBC_EUR_chr1 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr2 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr2.tsv.gz")
# sQTL_WBC_EUR_chr2 <- sQTL_WBC_EUR_chr2 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr2 <- sQTL_WBC_EUR_chr2 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr2 <- sQTL_WBC_EUR_chr2 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr2 <- sQTL_WBC_EUR_chr2 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr2 <- sQTL_WBC_EUR_chr2 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr3 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr3.tsv.gz")
# sQTL_WBC_EUR_chr3 <- sQTL_WBC_EUR_chr3 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr3 <- sQTL_WBC_EUR_chr3 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr3 <- sQTL_WBC_EUR_chr3 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr3 <- sQTL_WBC_EUR_chr3 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr3 <- sQTL_WBC_EUR_chr3 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr4 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr4.tsv.gz")
# sQTL_WBC_EUR_chr4 <- sQTL_WBC_EUR_chr4 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr4 <- sQTL_WBC_EUR_chr4 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr4 <- sQTL_WBC_EUR_chr4 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr4 <- sQTL_WBC_EUR_chr4 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr4 <- sQTL_WBC_EUR_chr4 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr5 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr5.tsv.gz")
# sQTL_WBC_EUR_chr5 <- sQTL_WBC_EUR_chr5 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr5 <- sQTL_WBC_EUR_chr5 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr5 <- sQTL_WBC_EUR_chr5 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr5 <- sQTL_WBC_EUR_chr5 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr5 <- sQTL_WBC_EUR_chr5 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr6 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr6.tsv.gz")
# sQTL_WBC_EUR_chr6 <- sQTL_WBC_EUR_chr6 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr6 <- sQTL_WBC_EUR_chr6 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr6 <- sQTL_WBC_EUR_chr6 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr6 <- sQTL_WBC_EUR_chr6 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr6 <- sQTL_WBC_EUR_chr6 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr7 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr7.tsv.gz")
# sQTL_WBC_EUR_chr7 <- sQTL_WBC_EUR_chr7 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr7 <- sQTL_WBC_EUR_chr7 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr7 <- sQTL_WBC_EUR_chr7 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr7 <- sQTL_WBC_EUR_chr7 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr7 <- sQTL_WBC_EUR_chr7 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr8 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr8.tsv.gz")
# sQTL_WBC_EUR_chr8 <- sQTL_WBC_EUR_chr8 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr8 <- sQTL_WBC_EUR_chr8 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr8 <- sQTL_WBC_EUR_chr8 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr8 <- sQTL_WBC_EUR_chr8 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr8 <- sQTL_WBC_EUR_chr8 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr9 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr9.tsv.gz")
# sQTL_WBC_EUR_chr9 <- sQTL_WBC_EUR_chr9 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr9 <- sQTL_WBC_EUR_chr9 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                           str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr9 <- sQTL_WBC_EUR_chr9 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr9 <- sQTL_WBC_EUR_chr9 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                     POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                     EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                     NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr9 <- sQTL_WBC_EUR_chr9 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr10 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr10.tsv.gz")
# sQTL_WBC_EUR_chr10 <- sQTL_WBC_EUR_chr10 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr10 <- sQTL_WBC_EUR_chr10 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr10 <- sQTL_WBC_EUR_chr10 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr10 <- sQTL_WBC_EUR_chr10 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr10 <- sQTL_WBC_EUR_chr10 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr11 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr11.tsv.gz")
# sQTL_WBC_EUR_chr11 <- sQTL_WBC_EUR_chr11 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr11 <- sQTL_WBC_EUR_chr11 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr11 <- sQTL_WBC_EUR_chr11 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr11 <- sQTL_WBC_EUR_chr11 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr11 <- sQTL_WBC_EUR_chr11 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr12 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr12.tsv.gz")
# sQTL_WBC_EUR_chr12 <- sQTL_WBC_EUR_chr12 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr12 <- sQTL_WBC_EUR_chr12 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr12 <- sQTL_WBC_EUR_chr12 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr12 <- sQTL_WBC_EUR_chr12 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr12 <- sQTL_WBC_EUR_chr12 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr13 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr13.tsv.gz")
# sQTL_WBC_EUR_chr13 <- sQTL_WBC_EUR_chr13 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr13 <- sQTL_WBC_EUR_chr13 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr13 <- sQTL_WBC_EUR_chr13 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr13 <- sQTL_WBC_EUR_chr13 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr13 <- sQTL_WBC_EUR_chr13 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr14 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr14.tsv.gz")
# sQTL_WBC_EUR_chr14 <- sQTL_WBC_EUR_chr14 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr14 <- sQTL_WBC_EUR_chr14 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr14 <- sQTL_WBC_EUR_chr14 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr14 <- sQTL_WBC_EUR_chr14 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr14 <- sQTL_WBC_EUR_chr14 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# # sQTL_WBC_EUR_chr14 %>% filter(variant_id == "chr14_24453482_G_T_b38" & grepl("chr14:24440926:24441674", phenotype_id))
#
# sQTL_WBC_EUR_chr15 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr15.tsv.gz")
# sQTL_WBC_EUR_chr15 <- sQTL_WBC_EUR_chr15 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr15 <- sQTL_WBC_EUR_chr15 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr15 <- sQTL_WBC_EUR_chr15 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr15 <- sQTL_WBC_EUR_chr15 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr15 <- sQTL_WBC_EUR_chr15 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr16 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr16.tsv.gz")
# sQTL_WBC_EUR_chr16 <- sQTL_WBC_EUR_chr16 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr16 <- sQTL_WBC_EUR_chr16 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr16 <- sQTL_WBC_EUR_chr16 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr16 <- sQTL_WBC_EUR_chr16 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr16 <- sQTL_WBC_EUR_chr16 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr17 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr17.tsv.gz")
# sQTL_WBC_EUR_chr17 <- sQTL_WBC_EUR_chr17 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr17 <- sQTL_WBC_EUR_chr17 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr17 <- sQTL_WBC_EUR_chr17 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr17 <- sQTL_WBC_EUR_chr17 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr17 <- sQTL_WBC_EUR_chr17 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr18 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr18.tsv.gz")
# sQTL_WBC_EUR_chr18 <- sQTL_WBC_EUR_chr18 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr18 <- sQTL_WBC_EUR_chr18 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr18 <- sQTL_WBC_EUR_chr18 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr18 <- sQTL_WBC_EUR_chr18 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr18 <- sQTL_WBC_EUR_chr18 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr19 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr19.tsv.gz")
# sQTL_WBC_EUR_chr19 <- sQTL_WBC_EUR_chr19 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr19 <- sQTL_WBC_EUR_chr19 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr19 <- sQTL_WBC_EUR_chr19 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr19 <- sQTL_WBC_EUR_chr19 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr19 <- sQTL_WBC_EUR_chr19 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr20 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr20.tsv.gz")
# sQTL_WBC_EUR_chr20 <- sQTL_WBC_EUR_chr20 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr20 <- sQTL_WBC_EUR_chr20 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr20 <- sQTL_WBC_EUR_chr20 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr20 <- sQTL_WBC_EUR_chr20 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr20 <- sQTL_WBC_EUR_chr20 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr21 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr21.tsv.gz")
# sQTL_WBC_EUR_chr21 <- sQTL_WBC_EUR_chr21 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr21 <- sQTL_WBC_EUR_chr21 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr21 <- sQTL_WBC_EUR_chr21 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr21 <- sQTL_WBC_EUR_chr21 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr21 <- sQTL_WBC_EUR_chr21 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR_chr22 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Whole_Blood.v8.EUR.sqtl_allpairs.chr22.tsv.gz")
# sQTL_WBC_EUR_chr22 <- sQTL_WBC_EUR_chr22 %>% filter(variant_id %in% sQTL_WBC$variant_id)
# sQTL_WBC_EUR_chr22 <- sQTL_WBC_EUR_chr22 %>% mutate(phenotype_id = paste0(str_split(phenotype_id, pattern="\\:", simplify = T)[,1],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,2],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,3],":",
#                                                                             str_split(phenotype_id, pattern="\\:", simplify = T)[,5]))
#
# sQTL_WBC_EUR_chr22 <- sQTL_WBC_EUR_chr22 %>% filter(phenotype_id %in% sQTL_WBC$phenotype_id)
#
# sQTL_WBC_EUR_chr22 <- sQTL_WBC_EUR_chr22 %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
#                                                       POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
#                                                       EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
#                                                       NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
#
# sQTL_WBC_EUR_chr22 <- sQTL_WBC_EUR_chr22 %>% inner_join(sQTL_WBC, by=c("phenotype_id"="phenotype_id", "variant_id"="variant_id"))
#
# sQTL_WBC_EUR <- bind_rows(sQTL_WBC_EUR_chr1, sQTL_WBC_EUR_chr2, sQTL_WBC_EUR_chr3, sQTL_WBC_EUR_chr4, sQTL_WBC_EUR_chr5,
#                            sQTL_WBC_EUR_chr6, sQTL_WBC_EUR_chr7, sQTL_WBC_EUR_chr8, sQTL_WBC_EUR_chr9, sQTL_WBC_EUR_chr10,
#                            sQTL_WBC_EUR_chr11, sQTL_WBC_EUR_chr12, sQTL_WBC_EUR_chr13, sQTL_WBC_EUR_chr14, sQTL_WBC_EUR_chr15,
#                            sQTL_WBC_EUR_chr16, sQTL_WBC_EUR_chr17, sQTL_WBC_EUR_chr18, sQTL_WBC_EUR_chr19, sQTL_WBC_EUR_chr20,
#                            sQTL_WBC_EUR_chr21, sQTL_WBC_EUR_chr22)

sQTL_WBC_EUR <- sQTL_WBC_EUR %>% mutate(N = round(ma_count.x/maf.x*(1/2)))
sQTL_WBC_EUR <- sQTL_WBC_EUR %>% mutate(SNP = paste0(CHR,":",POS))
exp_sQTL_WBC_EUR <- format_data(sQTL_WBC_EUR, type="exposure",
                                 phenotype_col = "phenotype_id",
                                 snp_col = "SNP",
                                 beta_col = "slope.x",
                                 se_col = "slope_se.x",
                                 effect_allele_col = "EA",
                                 other_allele_col = "NEA",
                                 eaf = "maf.x",
                                 pval_col = "pval_nominal.x",
                                 samplesize_col = "N",
                                 chr_col = "CHR",
                                 pos_col = "POS",
)
exp_sQTL_WBC_EUR <- exp_sQTL_WBC_EUR %>% filter(!(chr.exposure == "chr6" & pos.exposure > 28510120 & pos.exposure < 33480577))

# saveRDS(exp_sQTL_WBC_EUR, file="exposure_sQTL_WBC_EUR.rds")

exp_sQTL_WBC_EUR <- readRDS("exposure_sQTL_WBC_EUR.rds")
tmp <- readRDS("Whole_Blood_cluster_intron_gene_map.rds")
  
tmp <- tmp %>% mutate(exposure = paste0(pos,":",gene),
                      exposure1 = paste0(pos,":",cluster,":",gene)) %>% dplyr::select(exposure, exposure1)
exp_sQTL_WBC_EUR <- exp_sQTL_WBC_EUR %>% inner_join(tmp, by="exposure")
exp_sQTL_WBC_EUR <- exp_sQTL_WBC_EUR %>% mutate(exposure = exposure1)

saveRDS(exp_sQTL_WBC_EUR, file="exposure_sQTL_WBC_EUR.rds")


setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-eQTLMR/")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

eQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_Analysis_v8_eQTL_independent/Lung.v8.independent_eqtls.txt.gz")

eQTL_Lung_EUR_chr1 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr1.tsv.gz")
eQTL_Lung_EUR_chr1 <- eQTL_Lung_EUR_chr1 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr1 <- eQTL_Lung_EUR_chr1 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr1 <- eQTL_Lung_EUR_chr1 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr2 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr2.tsv.gz")
eQTL_Lung_EUR_chr2 <- eQTL_Lung_EUR_chr2 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr2 <- eQTL_Lung_EUR_chr2 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr2 <- eQTL_Lung_EUR_chr2 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr3 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr3.tsv.gz")
eQTL_Lung_EUR_chr3 <- eQTL_Lung_EUR_chr3 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr3 <- eQTL_Lung_EUR_chr3 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr3 <- eQTL_Lung_EUR_chr3 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr4 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr4.tsv.gz")
eQTL_Lung_EUR_chr4 <- eQTL_Lung_EUR_chr4 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr4 <- eQTL_Lung_EUR_chr4 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr4 <- eQTL_Lung_EUR_chr4 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr5 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr5.tsv.gz")
eQTL_Lung_EUR_chr5 <- eQTL_Lung_EUR_chr5 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr5 <- eQTL_Lung_EUR_chr5 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr5 <- eQTL_Lung_EUR_chr5 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr6 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr6.tsv.gz")
eQTL_Lung_EUR_chr6 <- eQTL_Lung_EUR_chr6 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr6 <- eQTL_Lung_EUR_chr6 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr6 <- eQTL_Lung_EUR_chr6 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr7 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr7.tsv.gz")
eQTL_Lung_EUR_chr7 <- eQTL_Lung_EUR_chr7 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr7 <- eQTL_Lung_EUR_chr7 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr7 <- eQTL_Lung_EUR_chr7 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr8 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr8.tsv.gz")
eQTL_Lung_EUR_chr8 <- eQTL_Lung_EUR_chr8 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr8 <- eQTL_Lung_EUR_chr8 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr8 <- eQTL_Lung_EUR_chr8 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr9 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr9.tsv.gz")
eQTL_Lung_EUR_chr9 <- eQTL_Lung_EUR_chr9 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr9 <- eQTL_Lung_EUR_chr9 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr9 <- eQTL_Lung_EUR_chr9 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr10 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr10.tsv.gz")
eQTL_Lung_EUR_chr10 <- eQTL_Lung_EUR_chr10 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr10 <- eQTL_Lung_EUR_chr10 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr10 <- eQTL_Lung_EUR_chr10 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr11 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr11.tsv.gz")
eQTL_Lung_EUR_chr11 <- eQTL_Lung_EUR_chr11 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr11 <- eQTL_Lung_EUR_chr11 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr11 <- eQTL_Lung_EUR_chr11 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr12 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr12.tsv.gz")
eQTL_Lung_EUR_chr12 <- eQTL_Lung_EUR_chr12 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr12 <- eQTL_Lung_EUR_chr12 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr12 <- eQTL_Lung_EUR_chr12 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr13 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr13.tsv.gz")
eQTL_Lung_EUR_chr13 <- eQTL_Lung_EUR_chr13 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr13 <- eQTL_Lung_EUR_chr13 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr13 <- eQTL_Lung_EUR_chr13 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr14 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr14.tsv.gz")
eQTL_Lung_EUR_chr14 <- eQTL_Lung_EUR_chr14 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr14 <- eQTL_Lung_EUR_chr14 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr14 <- eQTL_Lung_EUR_chr14 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr15 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr15.tsv.gz")
eQTL_Lung_EUR_chr15 <- eQTL_Lung_EUR_chr15 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr15 <- eQTL_Lung_EUR_chr15 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr15 <- eQTL_Lung_EUR_chr15 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr16 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr16.tsv.gz")
eQTL_Lung_EUR_chr16 <- eQTL_Lung_EUR_chr16 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr16 <- eQTL_Lung_EUR_chr16 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr16 <- eQTL_Lung_EUR_chr16 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr17 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr17.tsv.gz")
eQTL_Lung_EUR_chr17 <- eQTL_Lung_EUR_chr17 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr17 <- eQTL_Lung_EUR_chr17 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr17 <- eQTL_Lung_EUR_chr17 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr18 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr18.tsv.gz")
eQTL_Lung_EUR_chr18 <- eQTL_Lung_EUR_chr18 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr18 <- eQTL_Lung_EUR_chr18 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr18 <- eQTL_Lung_EUR_chr18 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr19 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr19.tsv.gz")
eQTL_Lung_EUR_chr19 <- eQTL_Lung_EUR_chr19 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr19 <- eQTL_Lung_EUR_chr19 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr19 <- eQTL_Lung_EUR_chr19 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr20 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr20.tsv.gz")
eQTL_Lung_EUR_chr20 <- eQTL_Lung_EUR_chr20 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr20 <- eQTL_Lung_EUR_chr20 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr20 <- eQTL_Lung_EUR_chr20 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr21 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr21.tsv.gz")
eQTL_Lung_EUR_chr21 <- eQTL_Lung_EUR_chr21 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr21 <- eQTL_Lung_EUR_chr21 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr21 <- eQTL_Lung_EUR_chr21 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Lung_EUR_chr22 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr22.tsv.gz")
eQTL_Lung_EUR_chr22 <- eQTL_Lung_EUR_chr22 %>% filter(variant_id %in% eQTL_Lung$variant_id)
eQTL_Lung_EUR_chr22 <- eQTL_Lung_EUR_chr22 %>% filter(phenotype_id %in% eQTL_Lung$gene_id)
eQTL_Lung_EUR_chr22 <- eQTL_Lung_EUR_chr22 %>% inner_join(eQTL_Lung, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))


eQTL_Lung_EUR <- bind_rows(eQTL_Lung_EUR_chr1, eQTL_Lung_EUR_chr2, eQTL_Lung_EUR_chr3, eQTL_Lung_EUR_chr4, eQTL_Lung_EUR_chr5,
                           eQTL_Lung_EUR_chr6, eQTL_Lung_EUR_chr7, eQTL_Lung_EUR_chr8, eQTL_Lung_EUR_chr9, eQTL_Lung_EUR_chr10,
                           eQTL_Lung_EUR_chr11, eQTL_Lung_EUR_chr12, eQTL_Lung_EUR_chr13, eQTL_Lung_EUR_chr14, eQTL_Lung_EUR_chr15,
                           eQTL_Lung_EUR_chr16, eQTL_Lung_EUR_chr17, eQTL_Lung_EUR_chr18, eQTL_Lung_EUR_chr19, eQTL_Lung_EUR_chr20,
                           eQTL_Lung_EUR_chr21, eQTL_Lung_EUR_chr22)

eQTL_Lung_EUR <- eQTL_Lung_EUR %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
                                                    POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                                    EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                                    NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])


eQTL_Lung_EUR <- eQTL_Lung_EUR %>% mutate(N = round(ma_count/maf.x*(1/2)))
eQTL_Lung_EUR <- eQTL_Lung_EUR %>% mutate(SNP = paste0(CHR,":",POS))
exp_eQTL_Lung_EUR <- format_data(eQTL_Lung_EUR, type="exposure",
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
exp_eQTL_Lung_EUR <- exp_eQTL_Lung_EUR %>% filter(!(chr.exposure == "chr6" & pos.exposure > 28510120 & pos.exposure < 33480577))

saveRDS(exp_eQTL_Lung_EUR, file="exposure_eQTL_Lung_EUR.rds")



eQTL_Whole_Blood <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_Analysis_v8_eQTL_independent/Whole_Blood.v8.independent_eqtls.txt.gz")

eQTL_Whole_Blood_EUR_chr1 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr1.tsv.gz")
eQTL_Whole_Blood_EUR_chr1 <- eQTL_Whole_Blood_EUR_chr1 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr1 <- eQTL_Whole_Blood_EUR_chr1 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr1 <- eQTL_Whole_Blood_EUR_chr1 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr2 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr2.tsv.gz")
eQTL_Whole_Blood_EUR_chr2 <- eQTL_Whole_Blood_EUR_chr2 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr2 <- eQTL_Whole_Blood_EUR_chr2 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr2 <- eQTL_Whole_Blood_EUR_chr2 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr3 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr3.tsv.gz")
eQTL_Whole_Blood_EUR_chr3 <- eQTL_Whole_Blood_EUR_chr3 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr3 <- eQTL_Whole_Blood_EUR_chr3 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr3 <- eQTL_Whole_Blood_EUR_chr3 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr4 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr4.tsv.gz")
eQTL_Whole_Blood_EUR_chr4 <- eQTL_Whole_Blood_EUR_chr4 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr4 <- eQTL_Whole_Blood_EUR_chr4 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr4 <- eQTL_Whole_Blood_EUR_chr4 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr5 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr5.tsv.gz")
eQTL_Whole_Blood_EUR_chr5 <- eQTL_Whole_Blood_EUR_chr5 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr5 <- eQTL_Whole_Blood_EUR_chr5 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr5 <- eQTL_Whole_Blood_EUR_chr5 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr6 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr6.tsv.gz")
eQTL_Whole_Blood_EUR_chr6 <- eQTL_Whole_Blood_EUR_chr6 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr6 <- eQTL_Whole_Blood_EUR_chr6 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr6 <- eQTL_Whole_Blood_EUR_chr6 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr7 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr7.tsv.gz")
eQTL_Whole_Blood_EUR_chr7 <- eQTL_Whole_Blood_EUR_chr7 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr7 <- eQTL_Whole_Blood_EUR_chr7 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr7 <- eQTL_Whole_Blood_EUR_chr7 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr8 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr8.tsv.gz")
eQTL_Whole_Blood_EUR_chr8 <- eQTL_Whole_Blood_EUR_chr8 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr8 <- eQTL_Whole_Blood_EUR_chr8 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr8 <- eQTL_Whole_Blood_EUR_chr8 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr9 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr9.tsv.gz")
eQTL_Whole_Blood_EUR_chr9 <- eQTL_Whole_Blood_EUR_chr9 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr9 <- eQTL_Whole_Blood_EUR_chr9 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr9 <- eQTL_Whole_Blood_EUR_chr9 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr10 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr10.tsv.gz")
eQTL_Whole_Blood_EUR_chr10 <- eQTL_Whole_Blood_EUR_chr10 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr10 <- eQTL_Whole_Blood_EUR_chr10 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr10 <- eQTL_Whole_Blood_EUR_chr10 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr11 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr11.tsv.gz")
eQTL_Whole_Blood_EUR_chr11 <- eQTL_Whole_Blood_EUR_chr11 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr11 <- eQTL_Whole_Blood_EUR_chr11 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr11 <- eQTL_Whole_Blood_EUR_chr11 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr12 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr12.tsv.gz")
eQTL_Whole_Blood_EUR_chr12 <- eQTL_Whole_Blood_EUR_chr12 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr12 <- eQTL_Whole_Blood_EUR_chr12 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr12 <- eQTL_Whole_Blood_EUR_chr12 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr13 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr13.tsv.gz")
eQTL_Whole_Blood_EUR_chr13 <- eQTL_Whole_Blood_EUR_chr13 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr13 <- eQTL_Whole_Blood_EUR_chr13 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr13 <- eQTL_Whole_Blood_EUR_chr13 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr14 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr14.tsv.gz")
eQTL_Whole_Blood_EUR_chr14 <- eQTL_Whole_Blood_EUR_chr14 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr14 <- eQTL_Whole_Blood_EUR_chr14 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr14 <- eQTL_Whole_Blood_EUR_chr14 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr15 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr15.tsv.gz")
eQTL_Whole_Blood_EUR_chr15 <- eQTL_Whole_Blood_EUR_chr15 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr15 <- eQTL_Whole_Blood_EUR_chr15 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr15 <- eQTL_Whole_Blood_EUR_chr15 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr16 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr16.tsv.gz")
eQTL_Whole_Blood_EUR_chr16 <- eQTL_Whole_Blood_EUR_chr16 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr16 <- eQTL_Whole_Blood_EUR_chr16 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr16 <- eQTL_Whole_Blood_EUR_chr16 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr17 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr17.tsv.gz")
eQTL_Whole_Blood_EUR_chr17 <- eQTL_Whole_Blood_EUR_chr17 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr17 <- eQTL_Whole_Blood_EUR_chr17 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr17 <- eQTL_Whole_Blood_EUR_chr17 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr18 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr18.tsv.gz")
eQTL_Whole_Blood_EUR_chr18 <- eQTL_Whole_Blood_EUR_chr18 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr18 <- eQTL_Whole_Blood_EUR_chr18 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr18 <- eQTL_Whole_Blood_EUR_chr18 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr19 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr19.tsv.gz")
eQTL_Whole_Blood_EUR_chr19 <- eQTL_Whole_Blood_EUR_chr19 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr19 <- eQTL_Whole_Blood_EUR_chr19 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr19 <- eQTL_Whole_Blood_EUR_chr19 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr20 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr20.tsv.gz")
eQTL_Whole_Blood_EUR_chr20 <- eQTL_Whole_Blood_EUR_chr20 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr20 <- eQTL_Whole_Blood_EUR_chr20 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr20 <- eQTL_Whole_Blood_EUR_chr20 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr21 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr21.tsv.gz")
eQTL_Whole_Blood_EUR_chr21 <- eQTL_Whole_Blood_EUR_chr21 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr21 <- eQTL_Whole_Blood_EUR_chr21 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr21 <- eQTL_Whole_Blood_EUR_chr21 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))

eQTL_Whole_Blood_EUR_chr22 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Whole_Blood.v8.EUR.allpairs.chr22.tsv.gz")
eQTL_Whole_Blood_EUR_chr22 <- eQTL_Whole_Blood_EUR_chr22 %>% filter(variant_id %in% eQTL_Whole_Blood$variant_id)
eQTL_Whole_Blood_EUR_chr22 <- eQTL_Whole_Blood_EUR_chr22 %>% filter(phenotype_id %in% eQTL_Whole_Blood$gene_id)
eQTL_Whole_Blood_EUR_chr22 <- eQTL_Whole_Blood_EUR_chr22 %>% inner_join(eQTL_Whole_Blood, by=c("phenotype_id"="gene_id", "variant_id"="variant_id"))


eQTL_Whole_Blood_EUR <- bind_rows(eQTL_Whole_Blood_EUR_chr1, eQTL_Whole_Blood_EUR_chr2, eQTL_Whole_Blood_EUR_chr3, eQTL_Whole_Blood_EUR_chr4, eQTL_Whole_Blood_EUR_chr5,
                           eQTL_Whole_Blood_EUR_chr6, eQTL_Whole_Blood_EUR_chr7, eQTL_Whole_Blood_EUR_chr8, eQTL_Whole_Blood_EUR_chr9, eQTL_Whole_Blood_EUR_chr10,
                           eQTL_Whole_Blood_EUR_chr11, eQTL_Whole_Blood_EUR_chr12, eQTL_Whole_Blood_EUR_chr13, eQTL_Whole_Blood_EUR_chr14, eQTL_Whole_Blood_EUR_chr15,
                           eQTL_Whole_Blood_EUR_chr16, eQTL_Whole_Blood_EUR_chr17, eQTL_Whole_Blood_EUR_chr18, eQTL_Whole_Blood_EUR_chr19, eQTL_Whole_Blood_EUR_chr20,
                           eQTL_Whole_Blood_EUR_chr21, eQTL_Whole_Blood_EUR_chr22)

eQTL_Whole_Blood_EUR <- eQTL_Whole_Blood_EUR %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1],
                                          POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                          EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                          NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])


eQTL_Whole_Blood_EUR <- eQTL_Whole_Blood_EUR %>% mutate(N = round(ma_count/maf.x*(1/2)))
eQTL_Whole_Blood_EUR <- eQTL_Whole_Blood_EUR %>% mutate(SNP = paste0(CHR,":",POS))
exp_eQTL_Whole_Blood_EUR <- format_data(eQTL_Whole_Blood_EUR, type="exposure",
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
exp_eQTL_Whole_Blood_EUR <- exp_eQTL_Whole_Blood_EUR %>% filter(!(chr.exposure == "chr6" & pos.exposure > 28510120 & pos.exposure < 33480577))

saveRDS(exp_eQTL_Whole_Blood_EUR, file="exposure_eQTL_WBC_EUR.rds")



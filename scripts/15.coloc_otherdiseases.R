setwd("~/scratch/09.COVID19/14.COVID19-esQTLMR/")

library(data.table)
library(coloc)
a2 <- fread("/project/richards/public/HGI/release7/pop_spec/COVID19_HGI_A2_ALL_eur_leave23andme_20220403_GRCh37.tsv.gz")
c2 <- fread("/project/richards/public/HGI/release7/pop_spec/COVID19_HGI_C2_ALL_eur_leave23andme_20220403_GRCh37.tsv.gz")

ipf <- fread("/project/richards/public/gwas/IPF/Allen_2019/allen_et_al_2019_ipf_meta_gwas_summary_statistics_public_download.txt")

#ATP11A

tmp <- a2 %>% filter(`#CHR` == 13 & POS > 113535741 - 500000 & POS < 113535741 + 500000)

tmp <- tmp %>% inner_join(ipf, by=c("#CHR"="chromosome", "POS"="position"))

tmp1 <- tmp %>% filter((b38_ref == non_effect_allele & b38_alt == effect_allele) | (b38_ref == effect_allele & b38_alt == non_effect_allele))

coloc.res <- coloc.abf(dataset1=list(beta=tmp1$all_inv_var_meta_beta, 
                                     varbeta = tmp1$all_inv_var_meta_sebeta*tmp1$all_inv_var_meta_sebeta,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF, 
                                     N=(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),
                                     s=tmp1$all_inv_var_meta_cases/(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(varbeta=tmp1$standard_error*tmp1$standard_error,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$effect_allele_frequency,
                                     beta=tmp1$beta,
                                     N=11259,
                                     s=0.236966,
                                     type="cc"))

# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 3.94e-17  1.73e-07  1.24e-12  4.47e-03  9.96e-01 
# [1] "PP abf for shared variant: 99.6%"

#DPP9

tmp <- a2 %>% filter(`#CHR` == 19 & POS > 4717672 - 500000 & POS < 4717672 + 500000)

tmp <- tmp %>% inner_join(ipf, by=c("#CHR"="chromosome", "POS"="position"))

tmp1 <- tmp %>% filter((b38_ref == non_effect_allele & b38_alt == effect_allele) | (b38_ref == effect_allele & b38_alt == non_effect_allele))


coloc.res <- coloc.abf(dataset1=list(beta=tmp1$all_inv_var_meta_beta, 
                                     varbeta = tmp1$all_inv_var_meta_sebeta*tmp1$all_inv_var_meta_sebeta,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF, 
                                     N=(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),
                                     s=tmp1$all_inv_var_meta_cases/(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(varbeta=tmp1$standard_error*tmp1$standard_error,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$effect_allele_frequency,
                                     beta=tmp1$beta,
                                     N=11259,
                                     s=0.236966,
                                     type="cc"))


# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 4.95e-53  3.32e-09  2.88e-47  9.34e-04  9.99e-01 
# [1] "PP abf for shared variant: 99.9%"


#NPNT
fev1fvc <- fread("/scratch/richards/tomoko.nakanishi/UKBBspirometry/GWAS/Shrine_30804560/Shrine_30804560_FEV1_to_FVC_RATIO_meta-analysis.txt.gz")

tmp <- a2 %>% filter(`#CHR` == 4 & POS > 106819053 - 500000 & POS < 106819053 + 500000)

tmp <- tmp %>% inner_join(fev1fvc, by=c("#CHR"="Chromosome", "POS"="Position_b37"))

tmp1 <- tmp %>% filter((b38_ref == Non_coded & b38_alt == Coded) | (b38_ref == Coded & b38_alt == Non_coded))

coloc.res <- coloc.abf(dataset1=list(beta=tmp1$all_inv_var_meta_beta, 
                                     varbeta = tmp1$all_inv_var_meta_sebeta*tmp1$all_inv_var_meta_sebeta,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF, 
                                     N=(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),
                                     s=tmp1$all_inv_var_meta_cases/(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(varbeta=tmp1$SE*tmp1$SE,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF,
                                     beta=tmp1$beta,
                                     N=321047,
                                     type="quant"))
# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 1.17e-130 8.56e-127  1.37e-07  4.27e-06  1.00e+00 
# [1] "PP abf for shared variant: 100%"

asthma <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GWAS/asthma2020/Asthma_Meta_results_191206.v2.txt.gz")
asthma <- asthma %>% mutate(Chromosome = as.numeric(gsub("chr", "", Chr)))

a2_b38 <- fread("/project/richards/public/HGI/release7/pop_spec/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz")

tmp <- a2_b38 %>% filter(`#CHR` == 4 & POS > 105897896 - 500000 & POS < 105897896 + 500000)

tmp <- tmp %>% inner_join(asthma, by=c("#CHR"="Chromosome", "POS"="Pos"))

tmp1 <- tmp %>% filter((REF == A0 & ALT == A1) | (REF == A1 & ALT == A0))
tmp1 <- tmp1 %>% mutate(Ntotal = case_when(`UKB-info` == "na" ~ 16247 + 346486,
                                      `IS-info` == "na" ~ 52942 + 355713,
                                      TRUE ~ 16247 + 346486 + 52942 + 355713),
                        Ncase = case_when(`UKB-info` == "na" ~ 16247,
                                          `IS-info` == "na" ~ 52942,
                                          TRUE ~ 16247 + 52942))

tmp1 <- tmp1 %>% filter(!is.na(all_meta_AF))
coloc.res <- coloc.abf(dataset1=list(beta=tmp1$all_inv_var_meta_beta, 
                                     varbeta = tmp1$all_inv_var_meta_sebeta*tmp1$all_inv_var_meta_sebeta,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF, 
                                     N=(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),
                                     s=tmp1$all_inv_var_meta_cases/(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF,
                                     beta=log(tmp1$`OR-A1`),
                                     pvalues = tmp1$P,
                                     N=(tmp1$Ntotal),
                                     s=tmp1$Ncase/tmp1$Ntotal,				  
                                     type="cc"))
# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 1.28e-09  9.42e-06  1.52e-07  1.19e-04  1.00e+00 
# [1] "PP abf for shared variant: 100%"


#MUC1

ibd <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GWAS/GCST004131/28067908-GCST004131-EFO_0003767-build37.f.tsv.gz")
ibd <- ibd %>% mutate(other_allele = toupper(other_allele),
                      effect_allele = toupper(effect_allele))

tmp <- c2 %>% filter(`#CHR` == 1 & POS > 155162067 - 500000 & POS < 155162067 + 500000)

tmp <- tmp %>% inner_join(ibd, by=c("#CHR"="chromosome", "POS"="base_pair_location"))

tmp1 <- tmp %>% filter((b38_ref == other_allele & b38_alt == effect_allele) | (b38_ref == effect_allele & b38_alt == other_allele))

coloc.res <- coloc.abf(dataset1=list(beta=tmp1$all_inv_var_meta_beta, 
                                     varbeta = tmp1$all_inv_var_meta_sebeta*tmp1$all_inv_var_meta_sebeta,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF, 
                                     N=(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),
                                     s=tmp1$all_inv_var_meta_cases/(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(varbeta=tmp1$standard_error*tmp1$standard_error,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF,
                                     beta=tmp1$beta,
                                     N=59957,
                                     s=25042/59957,
                                     type="cc"))

#gastric cancer

gc <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GWAS/hum0197.v3.EUR.GC.v1/GWASsummary_GaC_EUR_SakaueKanai2020.auto.txt.gz")

tmp <- c2 %>% filter(`#CHR` == 1 & POS > 155162067 - 500000 & POS < 155162067 + 500000)

tmp <- tmp %>% inner_join(gc, by=c("#CHR"="CHR", "POS"="POS"))

tmp1 <- tmp %>% filter((b38_ref == Allele1 & b38_alt == Allele2) | (b38_ref == Allele2 & b38_alt == Allele1))

coloc.res <- coloc.abf(dataset1=list(beta=tmp1$all_inv_var_meta_beta, 
                                     varbeta = tmp1$all_inv_var_meta_sebeta*tmp1$all_inv_var_meta_sebeta,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF, 
                                     N=(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),
                                     s=tmp1$all_inv_var_meta_cases/(tmp1$all_inv_var_meta_cases + tmp1$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(varbeta=tmp1$SE*tmp1$SE,
                                     snp=tmp1$SNP, 
                                     MAF=tmp1$all_meta_AF,
                                     beta=tmp1$BETA,
                                     N=476116,
                                     s=1029/476116,
                                     type="cc"))


#urate

urate <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GWAS/Urate_Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz", fill=TRUE)
map <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GWAS/Urate_Tin2019/b38.bed", header=F)

urate <- urate %>% inner_join(map, by=c("RSID"="V4"))

tmp <- c2 %>% filter(`#CHR` == 1 & POS > 155162067 - 500000 & POS < 155162067 + 500000)

tmp <- tmp %>% inner_join(urate, by=c("#CHR"="Chr", "POS"="V2"))

coloc.res <- coloc.abf(dataset1=list(beta=tmp$all_inv_var_meta_beta, 
                                     varbeta = tmp$all_inv_var_meta_sebeta*tmp$all_inv_var_meta_sebeta,
                                     snp=tmp$SNP, 
                                     MAF=tmp$all_meta_AF, 
                                     N=(tmp$all_inv_var_meta_cases + tmp$all_inv_var_meta_controls),
                                     s=tmp$all_inv_var_meta_cases/(tmp$all_inv_var_meta_cases + tmp$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(varbeta=tmp$StdErr*tmp$StdErr,
                                     snp=tmp$SNP, 
                                     MAF=tmp$all_meta_AF,
                                     beta=tmp$Effect,
                                     N=tmp$n_total_sum,
                                     type="quant"))

# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 2.90e-08  1.01e-03  2.87e-08  2.73e-08  9.99e-01 
# [1] "PP abf for shared variant: 99.9%"

## stroke
stroke <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/GWAS/Stroke2018/29531354-GCST006908-HP_0002140.h.tsv.gz")
c2 <- fread("/project/richards/public/HGI/release7/pop_spec/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz")
tmp <- c2 %>% filter(`#CHR` == 1 & POS > 156236330 - 500000 & POS < 156236330 + 500000)
tmp$`#CHR` <- as.character(tmp$`#CHR`)
tmp <- tmp %>% inner_join(stroke, by=c("#CHR"="hm_chrom", "POS"="hm_pos"))

coloc.res <- coloc.abf(dataset1=list(beta=tmp$all_inv_var_meta_beta, 
                                     varbeta = tmp$all_inv_var_meta_sebeta*tmp$all_inv_var_meta_sebeta,
                                     snp=tmp$SNP, 
                                     MAF=tmp$all_meta_AF, 
                                     N=(tmp$all_inv_var_meta_cases + tmp$all_inv_var_meta_controls),
                                     s=tmp$all_inv_var_meta_cases/(tmp$all_inv_var_meta_cases + tmp$all_inv_var_meta_controls),				  
                                     type="cc"), 
                       dataset2=list(varbeta=tmp$standard_error*tmp$standard_error,
                                     snp=tmp$SNP, 
                                     MAF=tmp$hm_effect_allele_frequency,
                                     beta=tmp$beta,
                                     N=440328,
                                     type="cc"))
# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 7.94e-05  3.00e-04  1.53e-01  5.78e-01  2.68e-01 
# [1] "PP abf for shared variant: 26.8%"


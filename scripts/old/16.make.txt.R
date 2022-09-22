setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

map <- fread("sQTL_WBC.map")
sig <- read.xlsx("WBC.highcoloc.xlsx")
out <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in c(1,3)){
  splice <- map$group_id[grepl(sig$ensembleID[i], map$group_id)]
  exp <- fread(paste0("sQTL_WBC/",splice,".tsv"))
  exp <- exp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                        POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                        EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                        NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
  exp <- exp %>% filter(POS > map$POS[grepl(sig$ensembleID[i], map$group_id)] - 500000 & POS < map$POS[grepl(sig$ensembleID[i], map$group_id)] + 500000)
  exp <- exp%>% mutate(CHRPOS = paste0(CHR,":",POS),
                       N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="EA", "REF"="NEA"))
  dat2 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="NEA", "REF"="EA"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":",simplify=TRUE)[,1], ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,2], ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,3],
                                          ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,5]))
  dat <- dat %>% filter(grepl(sig$junction[i], junction))
  
  DT_A2 <- dat %>% mutate(SNP = rsid,
                          mlog10p.o = -log10(all_inv_var_meta_p),
                          CHR = gsub("chr", "", CHR),
                          mbp = POS.y/1000000,
                          mlog10p.e = -log10(pval_nominal)
  ) %>% dplyr::select(SNP, mlog10p.o, CHR, mbp, mlog10p.e)
  ld <- fread(paste0("LD/",sig$rsid[i],".ld")) %>% dplyr::select(SNP_B, R2)
  DT_A2 <- DT_A2 %>% inner_join(ld, by=c("SNP"="SNP_B"))
  write.table(DT_A2, file=paste0("sQTL_WBC/DT_A2.",sig$ensembleID[i],".txt"), quote=F, col.names = T, row.names = F, sep="\t")
}

out <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in c(1,3)){
  splice <- map$group_id[grepl(sig$ensembleID[i], map$group_id)]
  exp <- fread(paste0("sQTL_WBC/",splice,".tsv"))
  exp <- exp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                        POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                        EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                        NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
  exp <- exp %>% filter(POS > map$POS[grepl(sig$ensembleID[i], map$group_id)] - 500000 & POS < map$POS[grepl(sig$ensembleID[i], map$group_id)] + 500000)
  exp <- exp%>% mutate(CHRPOS = paste0(CHR,":",POS),
                       N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="EA", "REF"="NEA"))
  dat2 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="NEA", "REF"="EA"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":",simplify=TRUE)[,1], ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,2], ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,3],
                                          ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,5]))
  dat <- dat %>% filter(grepl(sig$junction[i], junction))
  
  DT_A2 <- dat %>% mutate(SNP = rsid,
                          mlog10p.o = -log10(all_inv_var_meta_p),
                          CHR = gsub("chr", "", CHR),
                          mbp = POS.y/1000000,
                          mlog10p.e = -log10(pval_nominal)
                          ) %>% dplyr::select(SNP, mlog10p.o, CHR, mbp, mlog10p.e)
  ld <- fread(paste0("LD/",sig$rsid[i],".ld")) %>% dplyr::select(SNP_B, R2)
  DT_A2 <- DT_A2 %>% inner_join(ld, by=c("SNP"="SNP_B"))
  write.table(DT_A2, file=paste0("sQTL_WBC/DT_B2.",sig$ensembleID[i],".txt"), quote=F, col.names = T, row.names = F, sep="\t")
}

out <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in c(1,3)){
  splice <- map$group_id[grepl(sig$ensembleID[i], map$group_id)]
  exp <- fread(paste0("sQTL_WBC/",splice,".tsv"))
  exp <- exp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                        POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                        EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                        NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
  exp <- exp %>% filter(POS > map$POS[grepl(sig$ensembleID[i], map$group_id)] - 500000 & POS < map$POS[grepl(sig$ensembleID[i], map$group_id)] + 500000)
  exp <- exp%>% mutate(CHRPOS = paste0(CHR,":",POS),
                       N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="EA", "REF"="NEA"))
  dat2 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="NEA", "REF"="EA"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":",simplify=TRUE)[,1], ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,2], ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,3],
                                          ":", str_split(phenotype_id, pattern = ":",simplify=TRUE)[,5]))
  dat <- dat %>% filter(grepl(sig$junction[i], junction))
  
  DT_A2 <- dat %>% mutate(SNP = rsid,
                          mlog10p.o = -log10(all_inv_var_meta_p),
                          CHR = gsub("chr", "", CHR),
                          mbp = POS.y/1000000,
                          mlog10p.e = -log10(pval_nominal)
  ) %>% dplyr::select(SNP, mlog10p.o, CHR, mbp, mlog10p.e)
  ld <- fread(paste0("LD/",sig$rsid[i],".ld")) %>% dplyr::select(SNP_B, R2)
  DT_A2 <- DT_A2 %>% inner_join(ld, by=c("SNP"="SNP_B"))
  write.table(DT_A2, file=paste0("sQTL_WBC/DT_C2.",sig$ensembleID[i],".txt"), quote=F, col.names = T, row.names = F, sep="\t")
}


  

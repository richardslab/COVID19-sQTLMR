setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

map <- fread("sQTL_WBC.map")
out <- fread("/project/richards/public/HGI/release7/pop_spec/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in seq(1, dim(map)[1])){
  splice <- map$group_id[i]
  exp <- fread(paste0("sQTL_WBC/",splice,".tsv"))
  exp <- exp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                        POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                        EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                        NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
  exp <- exp %>% filter(POS > map$POS[i] - 500000 & POS < map$POS[i] + 500000)
  exp <- exp%>% mutate(CHRPOS = paste0(CHR,":",POS),
                       N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="EA", "REF"="NEA"))
  dat2 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="NEA", "REF"="EA"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(all_meta_AF = ifelse(is.na(all_meta_AF),maf, all_meta_AF))
  result <- data.frame(matrix(0, length(unique(dat$phenotype_id)), 7))
  for(j in seq(1, length(unique(dat$phenotype_id)))){
    datrev <- dat %>% filter(phenotype_id == unique(dat$phenotype_id)[j])
    coloc.res <- coloc.abf(dataset1=list(beta=datrev$slope, 
                                         varbeta = datrev$slope_se*datrev$slope_se,
                                         snp=datrev$variant_id, 
                                         MAF=datrev$maf, 
                                         N=datrev$N,
                                         type="quant"), 
                           dataset2=list(varbeta=datrev$all_inv_var_meta_sebeta*datrev$all_inv_var_meta_sebeta,
                                         snp=datrev$variant_id, 
                                         MAF=datrev$all_meta_AF,
                                         beta=datrev$all_inv_var_meta_beta,
                                         N=(datrev$all_inv_var_meta_cases + datrev$all_inv_var_meta_controls),
                                         s=datrev$all_inv_var_meta_cases/(datrev$all_inv_var_meta_cases + datrev$all_inv_var_meta_controls),					  
                                         type="cc"))
    result[j,1:7] <- c(unique(dat$phenotype_id)[j], coloc.res$summary)
  }
  write.table(result, file="A2.coloc.wbc.tsv", append=T, sep="\t", col.names = F, row.names = F, quote=F)
}

result <- fread("A2.coloc.wbc.tsv")

colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
write.table(result, file="A2.coloc.wbc.tsv", append=F, sep="\t", col.names = T, row.names = F, quote=F)


out <- fread("/project/richards/public/HGI/release7/pop_spec/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in seq(1, dim(map)[1])){
  splice <- map$group_id[i]
  exp <- fread(paste0("sQTL_WBC/",splice,".tsv"))
  exp <- exp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                        POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                        EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                        NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
  exp <- exp %>% filter(POS > map$POS[i] - 500000 & POS < map$POS[i] + 500000)
  exp <- exp%>% mutate(CHRPOS = paste0(CHR,":",POS),
                       N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="EA", "REF"="NEA"))
  dat2 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="NEA", "REF"="EA"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(all_meta_AF = ifelse(is.na(all_meta_AF),maf, all_meta_AF))
  result <- data.frame(matrix(0, length(unique(dat$phenotype_id)), 7))
  for(j in seq(1, length(unique(dat$phenotype_id)))){
    datrev <- dat %>% filter(phenotype_id == unique(dat$phenotype_id)[j])
    coloc.res <- coloc.abf(dataset1=list(beta=datrev$slope, 
                                         varbeta = datrev$slope_se*datrev$slope_se,
                                         snp=datrev$variant_id, 
                                         MAF=datrev$maf, 
                                         N=datrev$N,
                                         type="quant"), 
                           dataset2=list(varbeta=datrev$all_inv_var_meta_sebeta*datrev$all_inv_var_meta_sebeta,
                                         snp=datrev$variant_id, 
                                         MAF=datrev$all_meta_AF,
                                         beta=datrev$all_inv_var_meta_beta,
                                         N=(datrev$all_inv_var_meta_cases + datrev$all_inv_var_meta_controls),
                                         s=datrev$all_inv_var_meta_cases/(datrev$all_inv_var_meta_cases + datrev$all_inv_var_meta_controls),					  
                                         type="cc"))
    result[j,1:7] <- c(unique(dat$phenotype_id)[j], coloc.res$summary)
  }
  write.table(result, file="B2.coloc.wbc.tsv", append=T, sep="\t", col.names = F, row.names = F, quote=F)
}

result <- fread("B2.coloc.wbc.tsv")

colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
write.table(result, file="B2.coloc.wbc.tsv", append=F, sep="\t", col.names = T, row.names = F, quote=F)

out <- fread("/project/richards/public/HGI/release7/pop_spec/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in seq(1, dim(map)[1])){
  splice <- map$group_id[i]
  exp <- fread(paste0("sQTL_WBC/",splice,".tsv"))
  exp <- exp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                        POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                        EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                        NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
  exp <- exp %>% filter(POS > map$POS[i] - 500000 & POS < map$POS[i] + 500000)
  exp <- exp%>% mutate(CHRPOS = paste0(CHR,":",POS),
                       N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="EA", "REF"="NEA"))
  dat2 <- inner_join(out, exp, by=c("CHRPOS"="CHRPOS", "ALT"="NEA", "REF"="EA"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(all_meta_AF = ifelse(is.na(all_meta_AF),maf, all_meta_AF))
  result <- data.frame(matrix(0, length(unique(dat$phenotype_id)), 7))
  for(j in seq(1, length(unique(dat$phenotype_id)))){
    datrev <- dat %>% filter(phenotype_id == unique(dat$phenotype_id)[j])
    coloc.res <- coloc.abf(dataset1=list(beta=datrev$slope, 
                                         varbeta = datrev$slope_se*datrev$slope_se,
                                         snp=datrev$variant_id, 
                                         MAF=datrev$maf, 
                                         N=datrev$N,
                                         type="quant"), 
                           dataset2=list(varbeta=datrev$all_inv_var_meta_sebeta*datrev$all_inv_var_meta_sebeta,
                                         snp=datrev$variant_id, 
                                         MAF=datrev$all_meta_AF,
                                         beta=datrev$all_inv_var_meta_beta,
                                         N=(datrev$all_inv_var_meta_cases + datrev$all_inv_var_meta_controls),
                                         s=datrev$all_inv_var_meta_cases/(datrev$all_inv_var_meta_cases + datrev$all_inv_var_meta_controls),					  
                                         type="cc"))
    result[j,1:7] <- c(unique(dat$phenotype_id)[j], coloc.res$summary)
  }
  write.table(result, file="C2.coloc.wbc.tsv", append=T, sep="\t", col.names = F, row.names = F, quote=F)
}

result <- fread("C2.coloc.wbc.tsv")

colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
write.table(result, file="C2.coloc.wbc.tsv", append=F, sep="\t", col.names = T, row.names = F, quote=F)

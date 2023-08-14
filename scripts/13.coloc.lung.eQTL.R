setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

map <- fread("eQTL.map")
out <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

result <- data.frame(matrix(0, length(unique(map$exposure)), 7))
for(i in seq(1, dim(map)[1])){
  gene <- map$exposure[i]
  exp <- fread(paste0("eQTL_Lung/",gene,".tsv"))
  if(dim(exp)[1] >= 1){
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
    coloc.res <- coloc.abf(dataset1=list(beta=dat$slope, 
                                         varbeta = dat$slope_se*dat$slope_se,
                                         snp=dat$variant_id, 
                                         MAF=dat$maf, 
                                         N=dat$N,
                                         type="quant"), 
                           dataset2=list(varbeta=dat$all_inv_var_meta_sebeta*dat$all_inv_var_meta_sebeta,
                                         snp=dat$variant_id, 
                                         MAF=dat$all_meta_AF,
                                         beta=dat$all_inv_var_meta_beta,
                                         N=(dat$all_inv_var_meta_cases + dat$all_inv_var_meta_controls),
                                         s=dat$all_inv_var_meta_cases/(dat$all_inv_var_meta_cases + dat$all_inv_var_meta_controls),					  
                                         type="cc"))
    result[i,1:7] <- c(unique(dat$phenotype_id), coloc.res$summary) 
  }
}
colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

write.table(result, file="A2.coloc.Lung.eQTL.tsv", append=F, sep="\t", col.names = T, row.names = F, quote=F)

out <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in seq(1, dim(map)[1])){
  gene <- map$exposure[i]
  exp <- fread(paste0("eQTL_Lung/",gene,".tsv"))
  if(dim(exp)[1] >= 1){
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
    coloc.res <- coloc.abf(dataset1=list(beta=dat$slope, 
                                         varbeta = dat$slope_se*dat$slope_se,
                                         snp=dat$variant_id, 
                                         MAF=dat$maf, 
                                         N=dat$N,
                                         type="quant"), 
                           dataset2=list(varbeta=dat$all_inv_var_meta_sebeta*dat$all_inv_var_meta_sebeta,
                                         snp=dat$variant_id, 
                                         MAF=dat$all_meta_AF,
                                         beta=dat$all_inv_var_meta_beta,
                                         N=(dat$all_inv_var_meta_cases + dat$all_inv_var_meta_controls),
                                         s=dat$all_inv_var_meta_cases/(dat$all_inv_var_meta_cases + dat$all_inv_var_meta_controls),					  
                                         type="cc"))
    result[i,1:7] <- c(unique(dat$phenotype_id), coloc.res$summary) 
  }
}
colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
write.table(result, file="B2.coloc.Lung.eQTL.tsv", append=F, sep="\t", col.names = T, row.names = F, quote=F)


out <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/GWASsummary/pop_strat/COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz")
out <- out %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

for(i in seq(1, dim(map)[1])){
  gene <- map$exposure[i]
  exp <- fread(paste0("eQTL_Lung/",gene,".tsv"))
  if(dim(exp)[1] >= 1){
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
    coloc.res <- coloc.abf(dataset1=list(beta=dat$slope, 
                                         varbeta = dat$slope_se*dat$slope_se,
                                         snp=dat$variant_id, 
                                         MAF=dat$maf, 
                                         N=dat$N,
                                         type="quant"), 
                           dataset2=list(varbeta=dat$all_inv_var_meta_sebeta*dat$all_inv_var_meta_sebeta,
                                         snp=dat$variant_id, 
                                         MAF=dat$all_meta_AF,
                                         beta=dat$all_inv_var_meta_beta,
                                         N=(dat$all_inv_var_meta_cases + dat$all_inv_var_meta_controls),
                                         s=dat$all_inv_var_meta_cases/(dat$all_inv_var_meta_cases + dat$all_inv_var_meta_controls),					  
                                         type="cc"))
    result[i,1:7] <- c(unique(dat$phenotype_id), coloc.res$summary) 
  }
}
colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

write.table(result, file="C2.coloc.Lung.eQTL.tsv", append=F, sep="\t", col.names = T, row.names = F, quote=F)


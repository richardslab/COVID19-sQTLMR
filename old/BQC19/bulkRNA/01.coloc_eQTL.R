setwd("/home/tomoco/scratch/coloc/")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

args <- commandArgs(trailingOnly = TRUE)

covid <- fread(paste0("HGI/COVID19_HGI_",args[1],"_ALL_eur_leave23andme_20220403.tsv.gz"))

covid <- covid %>% mutate(CHRPOS = paste0("chr",`#CHR`,":",POS))

eQTL_noninf <- fread(paste0("eQTL/release6.4_2022-05-19_inf_No_pop_nfe.allpairs.chr",args[2],".awk.addedcolumns.txt"))
eQTL_inf <- fread(paste0("eQTL/release6.4_2022-05-19_inf_Yes_pop_nfe.allpairs.chr",args[2],".awk.addedcolumns.txt"))

eQTL_noninf1 <- eQTL_noninf %>% filter(POS > as.numeric(args[3]) - 500000 & POS < as.numeric(args[3]) + 500000)
eQTL_inf1 <- eQTL_inf %>% filter(POS > as.numeric(args[3]) - 500000 & POS < as.numeric(args[3]) + 500000)

result <- data.frame(matrix(0, length(unique(eQTL_noninf1$gene_id)), 7))

for(i in seq(1, length(unique(eQTL_noninf1$gene_id)))){
  expression <- unique(eQTL_noninf1$gene_id)[i]
  exp <- eQTL_noninf1 %>% filter(gene_id == expression)
  exp <- exp %>% mutate(CHRPOS = paste0(CHR,":",POS),
                        N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(covid, exp, by=c("CHRPOS"="CHRPOS", "ALT"="ALT", "REF"="REF"))
  dat2 <- inner_join(covid, exp, by=c("CHRPOS"="CHRPOS", "ALT"="REF", "REF"="ALT"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(all_meta_AF = ifelse(is.na(all_meta_AF),maf, all_meta_AF))
  datrev <- unique(dat) 
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
  result[i,1:7] <- c(expression, coloc.res$summary) 
}
colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

write.table(result, file=paste0("eQTL/",args[1],"_chr",args[2],"_",args[4],"_eQTL_noninfe.tsv"), sep="\t", col.names = T, row.names = F, quote=F)

result <- data.frame(matrix(0, length(unique(eQTL_inf1$gene_id)), 7))

for(i in seq(1, length(unique(eQTL_inf1$gene_id)))){
  expression <- unique(eQTL_inf1$gene_id)[i]
  exp <- eQTL_inf1 %>% filter(gene_id == expression)
  exp <- exp %>% mutate(CHRPOS = paste0(CHR,":",POS),
                        N = round(ma_count/maf*(1/2)))
  dat1 <- inner_join(covid, exp, by=c("CHRPOS"="CHRPOS", "ALT"="ALT", "REF"="REF"))
  dat2 <- inner_join(covid, exp, by=c("CHRPOS"="CHRPOS", "ALT"="REF", "REF"="ALT"))
  dat <- bind_rows(dat1, dat2)
  dat <- dat %>% mutate(all_meta_AF = ifelse(is.na(all_meta_AF),maf, all_meta_AF))
  datrev <- unique(dat) 
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
  result[i,1:7] <- c(expression, coloc.res$summary) 
}
colnames(result) <- c("phenotype_id", "NSNP","PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")

write.table(result, file=paste0("eQTL/",args[1],"_chr",args[2],"_",args[4],"_eQTL_infe.tsv"), sep="\t", col.names = T, row.names = F, quote=F)
            

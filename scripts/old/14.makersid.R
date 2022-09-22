setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

Lunglist <- read.xlsx("Lung.highcoloc.xlsx")
WBClist <- read.xlsx("WBC.highcoloc.xlsx")

rsid <- bind_rows(Lunglist, WBClist) %>% dplyr::select(rsid, CHR)
rsid <- unique(rsid)

write.table(rsid, file="sig_rsid", sep="\t", col.names = F, row.names = F, quote=F)

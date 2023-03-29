setwd("~/scratch/09.COVID19/14.COVID19-esQTLMR/LocusZoom")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr13.tsv.gz")
tmp <- tmp %>% dplyr::filter(phenotype_id == "ENSG00000068650.18")
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                  POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                  EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                  NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 112881427 - 500000 & POS < 112881427 + 500000)
write.table(tmp, file="ATP11A.eQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr4.tsv.gz")
tmp <- tmp %>% dplyr::filter(phenotype_id == "ENSG00000168743.12")
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 105897896 - 500000 & POS < 105897896 + 500000)
write.table(tmp, file="NPNT.eQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr19.tsv.gz")
tmp <- tmp %>% dplyr::filter(phenotype_id == "ENSG00000142002.16")
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 4717660 - 500000 & POS < 4717660 + 500000)
write.table(tmp, file="DPP9.eQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr1.tsv.gz")
tmp <- tmp %>% dplyr::filter(grepl("ENSG00000185499", phenotype_id))
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 155192276 - 500000 & POS < 155192276 + 500000)
write.table(tmp, file="MUC1.eQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr4.tsv.gz")
tmp <- tmp %>% dplyr::filter(phenotype_id == "ENSG00000168743.12")
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])



tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr13.tsv.gz")
tmp <- tmp %>% dplyr::filter(grepl("chr13:112875941:112880546", phenotype_id))
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 112881427 - 500000 & POS < 112881427 + 500000)
write.table(tmp, file="ATP11A.sQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr4.tsv.gz")
tmp <- tmp %>% dplyr::filter(grepl("chr4:105898001:105927336", phenotype_id))
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 105897896 - 500000 & POS < 105897896 + 500000)
write.table(tmp, file="NPNT.sQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr19.tsv.gz")
tmp <- tmp %>% dplyr::filter(grepl("chr19:4714337:4717615", phenotype_id))
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 4717660 - 500000 & POS < 4717660 + 500000)
write.table(tmp, file="DPP9.sQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")


tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx_summary/Lung.v8.EUR.allpairs.chr12.tsv.gz")
tmp <- tmp %>% dplyr::filter(phenotype_id == "ENSG00000089127.12")
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 112919388 - 500000 & POS < 112919388 + 500000)
write.table(tmp, file="OAS1.eQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")


tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr12.tsv.gz")
tmp <- tmp %>% dplyr::filter(grepl("chr12:112917700:112919389", phenotype_id))
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 112919388 - 500000 & POS < 112919388 + 500000)
write.table(tmp, file="OAS1.sQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr1.tsv.gz")
tmp <- tmp %>% dplyr::filter(grepl("chr1:155192310:155192786", phenotype_id))
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 155192276 - 500000 & POS < 155192276 + 500000)
write.table(tmp, file="MUC1.sQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")

tmp <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/summary/Lung.v8.EUR.sqtl_allpairs.chr1.tsv.gz")
tmp <- tmp %>% dplyr::filter(grepl("chr1:155205316:155206200", phenotype_id))
tmp <- tmp %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                      POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                      EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                      NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])

tmp <- tmp %>% dplyr::filter(POS > 155202588 - 500000 & POS < 155202588 + 500000)
write.table(tmp, file="THBS3.sQTL.tsv", quote=F, col.names = T, row.names = F, sep="\t")


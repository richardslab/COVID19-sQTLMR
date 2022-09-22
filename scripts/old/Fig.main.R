setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

lung <- read.xlsx("Lung.coloc.xlsx")
lunglist <- lung %>% filter(PP.H4.abf > PP.H0.abf + PP.H1.abf + PP.H2.abf + PP.H3.abf)
lunglist <- unique(lunglist$junction)

sQTL_Lung <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/sQTL/GTEx_Analysis_v8_sQTL_independent/Lung.v8.independent_sqtls.txt.gz")
sQTL_Lung <- sQTL_Lung %>% mutate(CHR = str_split(variant_id, pattern = "_",simplify=TRUE)[,1], 
                                  POS = as.numeric(str_split(variant_id, pattern = "_",simplify=TRUE)[,2]),
                                  EA = str_split(variant_id, pattern = "_",simplify=TRUE)[,4],
                                  NEA = str_split(variant_id, pattern = "_",simplify=TRUE)[,3])
sQTL_Lung <- sQTL_Lung %>% mutate(CHRPOS = paste0(CHR,":",POS),
                                  N = round(ma_count/maf*(1/2)))

#sQTL_Lung %>% filter(grepl("chr13:112875941:112880546:clu_3196:ENSG00000068650.18", phenotype_id))
mrresult <- read.xlsx("All_Lung_sQTL.xlsx")
mrresult_a <- mrresult %>% filter(outcome == "A2")#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_a <- mrresult_a %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_a <- unique(mrresult_a) %>% drop_na(p)
a2 <- fread("A2.coloc.lung.tsv")
a2 <- a2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
a2 <- unique(a2)
mrresult_a <- mrresult_a %>% left_join(a2, by="junction")
mrresult_a <- mrresult_a %>% filter(junction %in% lunglist)

mrresult_b <- mrresult %>% filter(outcome == "B2") #  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_b <- mrresult_b %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_b <- unique(mrresult_b) %>% drop_na(p)

b2 <- fread("B2.coloc.lung.tsv")
b2 <- b2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
mrresult_b <- mrresult_b %>% left_join(b2, by="junction")
mrresult_b <- mrresult_b %>% filter(junction %in% lunglist)

mrresult_c <- mrresult %>% filter(outcome == "C2")# %>% filter(p.adj < 0.05)#  %>% inner_join(a2, by=c("exposure"="phenotype_id"))
mrresult_c <- mrresult_c %>% mutate(junction = paste0(str_split(exposure, pattern = ":", simplify = T)[,1],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,2],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,3],":",
                                                      str_split(exposure, pattern = ":", simplify = T)[,5]))

mrresult_c <- unique(mrresult_c) %>% drop_na(p)
c2 <- fread("C2.coloc.lung.tsv")
c2 <- c2 %>% mutate(junction = paste0(str_split(phenotype_id, pattern = ":", simplify = T)[,1],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,2],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,3],":",
                                      str_split(phenotype_id, pattern = ":", simplify = T)[,5]))
mrresult_c <- mrresult_c %>% left_join(c2, by="junction")
mrresult_c <- mrresult_c %>% filter(junction %in% lunglist)

mrresult <- bind_rows(mrresult_a, mrresult_b, mrresult_c)

library(biomaRt)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(mrresult$gene))
colnames(gene) <- "ensembleID"

for(i in c(1:length(gene$ensembleID))){
  gene$gene_id[i] <- getBM(attributes='hgnc_symbol', filters = 'ensembl_gene_id',values = gene$ensembleID[i], mart = ensembl)
}

mrresult <- gene %>% right_join(mrresult, by=c("ensembleID"="gene"))

mrresult <- mrresult %>% mutate(splicing = paste0(gene_id," (", str_split(junction, pattern=":", simplify = T)[,1], ":",
                                                  str_split(junction, pattern=":", simplify = T)[,2], "-",
                                                  str_split(junction, pattern=":", simplify = T)[,3],")"))
                                                                                                              
tmp <- mrresult %>% filter(!grepl("chr1:", splicing))
tmp <- tmp %>% filter(!grepl("ABO", splicing))
library(ggplot2)
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/tmp.png",width=800, height=500)
tmp %>%
  ggplot(aes(y=OR, ymin=LL, ymax=UL, x=outcome, color=outcome)) + geom_pointrange(aes(color=outcome), lwd=0.8) + geom_hline(aes(), yintercept = 1, linetype=2) +
  theme_minimal() + facet_wrap(~ splicing, ncol = 1) + 
  scale_shape_identity() +  ggtitle("") +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=outcome), width=0.1, cex=1) + xlab("")  + ylab("Odds ratio (95%CI)") + 
  scale_color_brewer(palette = "Dark2") + 
  # scale_y_continuous(breaks=c(0.5, 1.0, 1.5, 2), limit=c(0.33, 2.0)) + 
  scale_fill_manual(values=c("#FF8000","#00994C")) + theme(axis.text.x = element_text(size=20,face="bold"),
                                                           axis.text.y=element_text(size=20,face="bold"),
                                                           axis.title=element_text(size=20,face="bold"), 
                                                           strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                           legend.title = element_blank(),
                                                           legend.text = element_blank(),
                                                           plot.title = element_text(size = 20),
                                                           strip.text = element_text(size=25)) +coord_flip() + guides(col = FALSE, fill = FALSE,
                                                                                                                      shape = FALSE)

dev.off()



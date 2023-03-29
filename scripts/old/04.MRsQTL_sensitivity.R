setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

All_Lung_sQTL <- readRDS("All_Lung_sQTL_ALL.rds")

tmp <- All_Lung_sQTL %>% filter(outcome == "Critical illness") %>% filter(p < 5e-8)

A2_Lung_sQTL_EUR <- readRDS("A2_Lung_sQTL_EUR.rds")
A2_Lung_sQTL_coloc <- fread("GTEx/sQTL/A2_sQTL_GTEx_Lung_coloc.tsv.gz")
A2_Lung_sQTL_EUR <- A2_Lung_sQTL_EUR %>% filter(exposure %in% tmp$exposure) %>% drop_na(p)
A2_Lung_sQTL_EUR <- A2_Lung_sQTL_EUR %>% left_join(A2_Lung_sQTL_coloc , by=c("exposure"="phenotype_id"))
A2_Lung_sQTL_EUR <- A2_Lung_sQTL_EUR %>% mutate(outcome = "Critical illness")

tmp <- All_Lung_sQTL %>% filter(outcome == "Hospitalization") %>% filter(p < 5e-8)

B2_Lung_sQTL_EUR <- readRDS("B2_Lung_sQTL_EUR.rds")
B2_Lung_sQTL_coloc <- fread("GTEx/sQTL/B2_sQTL_GTEx_Lung_coloc.tsv.gz")
B2_Lung_sQTL_EUR <- B2_Lung_sQTL_EUR %>% filter(exposure %in% tmp$exposure) %>% drop_na(p)
B2_Lung_sQTL_EUR <- B2_Lung_sQTL_EUR %>% left_join(B2_Lung_sQTL_coloc, by=c("exposure"="phenotype_id"))
B2_Lung_sQTL_EUR <- B2_Lung_sQTL_EUR %>% mutate(outcome = "Hospitalization")

tmp <- All_Lung_sQTL %>% filter(outcome == "Reported infection") %>% filter(p < 5e-8)

C2_Lung_sQTL_EUR <- readRDS("C2_Lung_sQTL_EUR.rds")
C2_Lung_sQTL_coloc <- fread("GTEx/sQTL/C2_sQTL_GTEx_Lung_coloc.tsv.gz")
C2_Lung_sQTL_EUR <- C2_Lung_sQTL_EUR %>% filter(exposure %in% tmp$exposure) %>% drop_na(p)
C2_Lung_sQTL_EUR <- C2_Lung_sQTL_EUR %>% left_join(C2_Lung_sQTL_coloc, by=c("exposure"="phenotype_id"))
C2_Lung_sQTL_EUR <- C2_Lung_sQTL_EUR %>% mutate(outcome = "Reported infection")

ALL_Lung_sQTL_EUR <- bind_rows(A2_Lung_sQTL_EUR, B2_Lung_sQTL_EUR, C2_Lung_sQTL_EUR)

ALL_Lung_sQTL_EUR <- ALL_Lung_sQTL_EUR %>% mutate(OR = exp(b),
         LL = exp(b + qnorm(0.025)*se),
         UL = exp(b + qnorm(0.975)*se),
         ensembleID = str_split(exposure, pattern="\\:", simplify = T)[,4],
         ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(ALL_Lung_sQTL_EUR$ensembleID))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

ALL_Lung_sQTL_EUR <- ALL_Lung_sQTL_EUR %>% left_join(G_list, by=c("ensembleID"="ensembl_gene_id"))

ALL_Lung_sQTL_EUR %>% saveRDS("ALL_Lung_sQTL_EUR_coloc.rds")

final <- ALL_Lung_sQTL_EUR %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, SNP, OR, LL, UL, p,
                                             NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

final %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable3.Lung.sensitivity.xlsx")

ALL_Lung_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
ALL_Lung_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
ALL_Lung_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()

All_WBC_sQTL <- readRDS("All_WBC_sQTL_ALL.rds")

tmp <- All_WBC_sQTL %>% filter(outcome == "Critical illness") %>% filter(p < 5e-8)

A2_WBC_sQTL_EUR <- readRDS("A2_WBC_sQTL_EUR.rds")
A2_WBC_sQTL_coloc <- fread("GTEx/sQTL/A2_sQTL_GTEx_Whole_Blood_coloc.tsv.gz")
A2_WBC_sQTL_EUR <- A2_WBC_sQTL_EUR %>% filter(exposure %in% tmp$exposure) %>% drop_na(p)
A2_WBC_sQTL_EUR <- A2_WBC_sQTL_EUR %>% left_join(A2_WBC_sQTL_coloc , by=c("exposure"="phenotype_id"))
A2_WBC_sQTL_EUR <- A2_WBC_sQTL_EUR %>% mutate(outcome = "Critical illness")

tmp <- All_WBC_sQTL %>% filter(outcome == "Hospitalization") %>% filter(p < 5e-8)

B2_WBC_sQTL_EUR <- readRDS("B2_WBC_sQTL_EUR.rds")
B2_WBC_sQTL_coloc <- fread("GTEx/sQTL/B2_sQTL_GTEx_Whole_Blood_coloc.tsv.gz")
B2_WBC_sQTL_EUR <- B2_WBC_sQTL_EUR %>% filter(exposure %in% tmp$exposure) %>% drop_na(p)
B2_WBC_sQTL_EUR <- B2_WBC_sQTL_EUR %>% left_join(B2_WBC_sQTL_coloc, by=c("exposure"="phenotype_id"))
B2_WBC_sQTL_EUR <- B2_WBC_sQTL_EUR %>% mutate(outcome = "Hospitalization")

tmp <- All_WBC_sQTL %>% filter(outcome == "Reported infection") %>% filter(p < 5e-8)

C2_WBC_sQTL_EUR <- readRDS("C2_WBC_sQTL_EUR.rds")
C2_WBC_sQTL_coloc <- fread("GTEx/sQTL/C2_sQTL_GTEx_Whole_Blood_coloc.tsv.gz")
C2_WBC_sQTL_EUR <- C2_WBC_sQTL_EUR %>% filter(exposure %in% tmp$exposure) %>% drop_na(p)
C2_WBC_sQTL_EUR <- C2_WBC_sQTL_EUR %>% left_join(C2_WBC_sQTL_coloc, by=c("exposure"="phenotype_id"))
C2_WBC_sQTL_EUR <- C2_WBC_sQTL_EUR %>% mutate(outcome = "Reported infection")

ALL_WBC_sQTL_EUR <- bind_rows(A2_WBC_sQTL_EUR, B2_WBC_sQTL_EUR, C2_WBC_sQTL_EUR)

ALL_WBC_sQTL_EUR <- ALL_WBC_sQTL_EUR %>% mutate(OR = exp(b),
                                                  LL = exp(b + qnorm(0.025)*se),
                                                  UL = exp(b + qnorm(0.975)*se),
                                                  ensembleID = str_split(exposure, pattern="\\:", simplify = T)[,4],
                                                  ensembleID = str_split(ensembleID, pattern="\\.", simplify = T)[,1])

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
gene <- data.frame(unique(ALL_WBC_sQTL_EUR$ensembleID))
colnames(gene) <- "ensembleID"

G_list <- getBM(attributes=c("ensembl_gene_id", 'hgnc_symbol'), filters = 'ensembl_gene_id',values = gene$ensembleID,mart = ensembl)

ALL_WBC_sQTL_EUR <- ALL_WBC_sQTL_EUR %>% left_join(G_list, by=c("ensembleID"="ensembl_gene_id"))

ALL_WBC_sQTL_EUR %>% saveRDS("ALL_WBC_sQTL_EUR_coloc.rds")

final <- ALL_WBC_sQTL_EUR %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, SNP, OR, LL, UL, p,
                                             NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

final %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable4.WBC.sensitivity.xlsx")

ALL_WBC_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
ALL_WBC_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
ALL_WBC_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()

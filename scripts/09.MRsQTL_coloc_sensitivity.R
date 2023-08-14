setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

tmp <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable1.Lung.MR.xlsx") %>% 
filter(p < 0.05/27230)
tmp <- unique(tmp$exposure)

A2_Lung_sQTL_EUR <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable1.Lung.MR.xlsx") 
A2_Lung_sQTL_EUR <- A2_Lung_sQTL_EUR %>% filter(outcome == "Critical illness") %>% filter(exposure %in% tmp)
#length(unique(A2_Lung_sQTL_EUR$exposure)) 7
A2_Lung_sQTL_coloc <- fread("GTEx/sQTL/A2_sQTL_GTEx_Lung_coloc.tsv.gz")

A2_Lung_sQTL_EUR <- A2_Lung_sQTL_EUR %>% mutate(phenotype_id = paste0(str_split(exposure, pattern=":", simplify = T)[,1],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,2],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,3],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,5]))
A2_Lung_sQTL_EUR <- A2_Lung_sQTL_EUR %>% drop_na(p)
A2_Lung_sQTL_EUR <- A2_Lung_sQTL_EUR %>% left_join(A2_Lung_sQTL_coloc , by=c("phenotype_id"))

B2_Lung_sQTL_EUR <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable1.Lung.MR.xlsx") 
B2_Lung_sQTL_EUR <- B2_Lung_sQTL_EUR %>% filter(outcome == "Hospitalization") %>% filter(exposure %in% tmp)
#length(unique(B2_Lung_sQTL_EUR$exposure)) 7
B2_Lung_sQTL_coloc <- fread("GTEx/sQTL/B2_sQTL_GTEx_Lung_coloc.tsv.gz")

B2_Lung_sQTL_EUR <- B2_Lung_sQTL_EUR %>% mutate(phenotype_id = paste0(str_split(exposure, pattern=":", simplify = T)[,1],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,2],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,3],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,5]))
B2_Lung_sQTL_EUR <- B2_Lung_sQTL_EUR %>% drop_na(p)
B2_Lung_sQTL_EUR <- B2_Lung_sQTL_EUR %>% left_join(B2_Lung_sQTL_coloc , by=c("phenotype_id"))


C2_Lung_sQTL_EUR <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable1.Lung.MR.xlsx") 
C2_Lung_sQTL_EUR <- C2_Lung_sQTL_EUR %>% filter(outcome == "Reported infection") %>% filter(exposure %in% tmp)
#length(unique(C2_Lung_sQTL_EUR$exposure)) 7
C2_Lung_sQTL_coloc <- fread("GTEx/sQTL/C2_sQTL_GTEx_Lung_coloc.tsv.gz")

C2_Lung_sQTL_EUR <- C2_Lung_sQTL_EUR %>% mutate(phenotype_id = paste0(str_split(exposure, pattern=":", simplify = T)[,1],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,2],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,3],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,5]))
C2_Lung_sQTL_EUR <- C2_Lung_sQTL_EUR %>% drop_na(p)
C2_Lung_sQTL_EUR <- C2_Lung_sQTL_EUR %>% left_join(C2_Lung_sQTL_coloc , by=c("phenotype_id"))

EUR_Lung_sQTL_EUR <- bind_rows(A2_Lung_sQTL_EUR, B2_Lung_sQTL_EUR, C2_Lung_sQTL_EUR)

EUR_Lung_sQTL_EUR %>% saveRDS("EUR_Lung_sQTL_EUR_coloc.rds")

final <- EUR_Lung_sQTL_EUR %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, 
                                             NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

final %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable4.Lung.sensitivity.xlsx")

EUR_Lung_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()#5
EUR_Lung_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()#4
EUR_Lung_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()#8


tmp <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable2.WBC.MR.xlsx") %>% 
  filter(p < 0.05/27230)
tmp <- unique(tmp$exposure)

A2_WBC_sQTL_EUR <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable2.WBC.MR.xlsx") 
A2_WBC_sQTL_EUR <- A2_WBC_sQTL_EUR %>% filter(outcome == "Critical illness") %>% filter(exposure %in% tmp)
#length(unique(A2_WBC_sQTL_EUR$exposure)) 7
A2_WBC_sQTL_coloc <- fread("GTEx/sQTL/A2_sQTL_GTEx_Whole_Blood_coloc.tsv.gz")

A2_WBC_sQTL_EUR <- A2_WBC_sQTL_EUR %>% mutate(phenotype_id = paste0(str_split(exposure, pattern=":", simplify = T)[,1],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,2],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,3],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,5]))
A2_WBC_sQTL_EUR <- A2_WBC_sQTL_EUR %>% drop_na(p)
A2_WBC_sQTL_EUR <- A2_WBC_sQTL_EUR %>% left_join(A2_WBC_sQTL_coloc , by=c("phenotype_id"))

B2_WBC_sQTL_EUR <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable2.WBC.MR.xlsx") 
B2_WBC_sQTL_EUR <- B2_WBC_sQTL_EUR %>% filter(outcome == "Hospitalization") %>% filter(exposure %in% tmp)
#length(unique(B2_WBC_sQTL_EUR$exposure)) 7
B2_WBC_sQTL_coloc <- fread("GTEx/sQTL/B2_sQTL_GTEx_Whole_Blood_coloc.tsv.gz")

B2_WBC_sQTL_EUR <- B2_WBC_sQTL_EUR %>% mutate(phenotype_id = paste0(str_split(exposure, pattern=":", simplify = T)[,1],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,2],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,3],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,5]))
B2_WBC_sQTL_EUR <- B2_WBC_sQTL_EUR %>% drop_na(p)
B2_WBC_sQTL_EUR <- B2_WBC_sQTL_EUR %>% left_join(B2_WBC_sQTL_coloc , by=c("phenotype_id"))


C2_WBC_sQTL_EUR <- readxl::read_excel("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable2.WBC.MR.xlsx") 
C2_WBC_sQTL_EUR <- C2_WBC_sQTL_EUR %>% filter(outcome == "Reported infection") %>% filter(exposure %in% tmp)
#length(unique(C2_WBC_sQTL_EUR$exposure)) 7
C2_WBC_sQTL_coloc <- fread("GTEx/sQTL/C2_sQTL_GTEx_Whole_Blood_coloc.tsv.gz")

C2_WBC_sQTL_EUR <- C2_WBC_sQTL_EUR %>% mutate(phenotype_id = paste0(str_split(exposure, pattern=":", simplify = T)[,1],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,2],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,3],":",
                                                                      str_split(exposure, pattern=":", simplify = T)[,5]))
C2_WBC_sQTL_EUR <- C2_WBC_sQTL_EUR %>% drop_na(p)
C2_WBC_sQTL_EUR <- C2_WBC_sQTL_EUR %>% left_join(C2_WBC_sQTL_coloc , by=c("phenotype_id"))

EUR_WBC_sQTL_EUR <- bind_rows(A2_WBC_sQTL_EUR, B2_WBC_sQTL_EUR, C2_WBC_sQTL_EUR)

EUR_WBC_sQTL_EUR %>% saveRDS("EUR_WBC_sQTL_EUR_coloc.rds")

final <- EUR_WBC_sQTL_EUR %>% dplyr::select(exposure, ensembleID, hgnc_symbol, outcome, 
                                             NSNP, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)

final %>% write.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable5.WBC.sensitivity.xlsx")

EUR_WBC_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8) %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()#1
EUR_WBC_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8) %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()#1
EUR_WBC_sQTL_EUR %>% filter(p < 0.05/27230 & PP.H4.abf > 0.8) %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()#3

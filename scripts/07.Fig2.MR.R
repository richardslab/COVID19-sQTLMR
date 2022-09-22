setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

ALL_Lung_sQTL_EUR <- readRDS("ALL_Lung_sQTL_EUR_coloc.rds")
ALL_WBC_sQTL_EUR <- readRDS("ALL_WBC_sQTL_EUR_coloc.rds")

tmp_lung <- ALL_Lung_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8)
tmp_lung %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_lung %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_lung %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()

tmp_wbc <- ALL_WBC_sQTL_EUR %>% filter(p < 0.05 & PP.H4.abf > 0.8)
tmp_wbc %>% filter(outcome == "Critical illness") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_wbc %>% filter(outcome == "Hospitalization") %>% dplyr::select(exposure) %>% unique() %>% dim()
tmp_wbc %>% filter(outcome == "Reported infection") %>% dplyr::select(exposure) %>% unique() %>% dim()


tmp <- bind_rows(tmp_lung, tmp_wbc)

Lung <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable1.Lung.MR.xlsx")
Lung <- Lung %>% filter(exposure %in% tmp_lung$exposure)

expression_Lung <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable5.eQTL_Lung.xlsx")
expression_Lung <- expression_Lung %>% filter(PP.H4.abf > 0.8)

expression_WBC <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable6.eQTL_WBC.xlsx")
expression_WBC <- expression_WBC %>% filter(PP.H4.abf.GTEx > 0.8)

expression <- unique(c(expression_WBC$gene, expression_Lung$gene))

Lung <- Lung %>% filter(hgnc_symbol != "")  %>% filter(!(ensembleID %in% expression))
Lung <- Lung %>% filter(!(hgnc_symbol %in% c( "THBS3")))

Lung <- Lung %>% mutate(splicing = paste0(hgnc_symbol," (", str_split(exposure, pattern=":", simplify = T)[,1], ":",
                                                  str_split(exposure, pattern=":", simplify = T)[,2], "-",
                                                  str_split(exposure, pattern=":", simplify = T)[,3],")"))

library(ggplot2)
library(dplyr)
library(forcats)
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/Fig2.png",width=800, height=800)
Lung %>% mutate(outcome = fct_relevel(outcome, 
                                      "Reported infection", "Hospitalization", "Critical illness")) %>%
  mutate(splicing = fct_relevel(splicing, 
                                "ATP11A (chr13:112875941-112880546)","DPP9 (chr19:4714337-4717615)",  "NPNT (chr4:105898001-105927336)","OAS1 (chr12:112917700-112919389)",
                                "MUC1 (chr1:155192310-155192786)")) %>%
  ggplot(aes(y=OR, ymin=LL, ymax=UL, x=outcome, color=outcome)) + geom_pointrange(aes(color=outcome), lwd=0.8) + geom_hline(aes(), yintercept = 1, linetype=2) +
  theme_minimal() + facet_wrap(~ splicing, ncol = 1) + 
  scale_shape_identity() +  ggtitle("") +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=outcome), width=0.1, cex=1) + xlab("")  + ylab("Odds ratio (95%CI) per normalised read counts") + 
  scale_color_manual(values=c("#00E6E6", "#FFD21A", "#FF9A00")) + 
  #scale_y_continuous(breaks=c(0.5, 1.0, 1.5, 2), limit=c(0.33, 2.0)) + 
  scale_fill_manual(values=c("#00E6E6", "#FFD21A", "#FF9A00")) + theme(axis.text.x = element_text(size=20,face="bold"),
                                                           axis.text.y=element_text(size=20,face="bold"),
                                                           axis.title=element_text(size=20,face="bold"), 
                                                           strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=20,face="bold"),
                                                           legend.title = element_blank(),
                                                           legend.text = element_blank(),
                                                           plot.title = element_text(size = 20),
                                                           strip.text = element_text(size=25)) +coord_flip() + guides(col = FALSE, fill = FALSE,
                                                                                                                      shape = FALSE)
dev.off()


Whole_Blood <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/SuppleTable2.WBC.MR.xlsx")
Whole_Blood <- Whole_Blood %>% filter(exposure %in% tmp_wbc$exposure)


Whole_Blood <- Whole_Blood %>% filter(!(ensembleID %in% expression)) %>% filter(hgnc_symbol != "")

Whole_Blood <- Whole_Blood %>% mutate(splicing = paste0(hgnc_symbol," (", str_split(exposure, pattern=":", simplify = T)[,1], ":",
                                          str_split(exposure, pattern=":", simplify = T)[,2], "-",
                                          str_split(exposure, pattern=":", simplify = T)[,3],")"))

png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/tmp.png",width=800, height=1500)
Whole_Blood %>%
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


#ATG16L2 only critical both lung and whole blood
#BTN2A1 only hospital both lung and whole blood
#OAS1 severity + susceptibility both lung and whole blood

#ATP11A, DPP9, RAPH1  lung specific severity
#CEP85, WBC only severity

#EPS8L1, SLC11A1, lung specific only susceptibility

#GBAP1, THBS3, MUC1, PMF1

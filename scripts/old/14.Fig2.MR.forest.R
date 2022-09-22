setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR")

library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(coloc)

mrresult <- read.xlsx("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/supTable3.coloc.lung.xlsx")
mrresult <- mrresult %>% mutate(splicing = paste0(gene_id," (", str_split(exposure, pattern=":", simplify = T)[,1], ":",
                                                  str_split(exposure, pattern=":", simplify = T)[,2], "-",
                                                  str_split(exposure, pattern=":", simplify = T)[,3],")"))
                                                      
tmp <- mrresult %>% filter(gene_id %in% c("ATP11A", "ABO", "NPNT",
                                          "DPP9", "SFTPA2", "OAS1"))                                                        
tmp <- tmp %>% filter(!grepl("chr9:133256356-133257409", splicing))
tmp <- tmp %>% filter(!grepl("chr19:4714337-4719851", splicing))
library(ggplot2)
library(forcats)
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/Fig2.MR.png",width=1000, height=800)
tmp %>%
  mutate(covid19_outcome = fct_relevel(outcome, 
                            "reported SARS-CoV-2 infection", "hospitalization", "critical illness")) %>%
  ggplot(aes(y=OR, ymin=LL, ymax=UL, x=covid19_outcome, color=covid19_outcome)) + geom_pointrange(aes(color=covid19_outcome), lwd=0.8) + geom_hline(aes(), yintercept = 1, linetype=2) +
  theme_minimal() + facet_wrap(~ splicing, ncol = 1) + 
  scale_shape_identity() +  ggtitle("") +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=covid19_outcome), width=0.1, cex=1) + xlab("")  + ylab("Odds ratio (95%CI)") + 
  scale_color_brewer(palette = "Set2") + 
  scale_y_continuous(breaks=c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3), limit=c(0.13, 1.3)) + 
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



setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/HPA/")

rna_single_cell_type_tissue <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/HPA/rna_single_cell_type_tissue.tsv.gz")

all <- rna_single_cell_type_tissue %>% filter(Tissue %in% c("lung", "pbmc")) %>% filter(`Gene name` %in% c("ATP11A", "DPP9", "NPNT", "MUC1", "OAS1"))
all <- all %>% mutate(celltype = paste0(`Cell type`," (", Tissue,": ",Cluster,")"))
all <- all %>% mutate(logTPM = log10(pTPM+1))

library("corrplot")

all <- all %>% arrange(Tissue)
tmp <- all[,c("Gene name", "celltype","logTPM", "Tissue")]
tmp <- tmp %>% arrange(Tissue, desc(logTPM))
tmp <- tmp %>% dplyr::select(-Tissue)

tmp1 <- reshape(tmp, idvar = "Gene name", timevar = "celltype", direction = "wide")


M <- as.matrix(tmp1[,-1])
rownames(M) <- tmp1$`Gene name`
colnames(M) <- str_split(colnames(tmp1)[-1], pattern="\\.", simplify = T)[,2]

png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/Fig4B.png",width=800, height=500)
corrplot(M, is.corr= FALSE, col=COL1(sequential = c("Red"), n=200), tl.col="black")
dev.off()
normal_tissue <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/HPA/normal_tissue.tsv.gz")
all <- normal_tissue %>% filter(`Gene name` %in% c("ATP11A", "DPP9", "NPNT", "MUC1", "OAS1"))
all <- all %>% mutate(TissueCelltype = paste0(Tissue," (", `Cell type`,")"))
high_medium_tissue <- unique(all$TissueCelltype[all$Level %in% c("High")])
all <- all %>% filter(TissueCelltype %in% high_medium_tissue)

library("corrplot")

all <- all %>% arrange(Tissue)
tmp <- all[,c("Gene name", "TissueCelltype","Level", "Tissue", "Cell type")]
library(forcats)
tmp <- tmp %>% mutate(level = fct_relevel(Level, 
                          "High", "Medium", "Low", 
                          "Not detected"),
                      tissue = fct_relevel(Tissue, 
                                           "lung", unique(all$Tissue)[-22]),
                      genename = fct_relevel(`Gene name`, 
                                             "ATP11A","DPP9", "NPNT", "OAS1", "MUC1"))
tmp <- tmp %>% arrange(tissue, level)

tmp <- tmp %>% mutate(tissuecelltype = fct_relevel(TissueCelltype, 
                                                   rev(unique(tmp$TissueCelltype))))

png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/Fig4A.png",width=400, height=500)
ggplot(tmp, aes(x=genename, y=tissuecelltype, fill=level)) + geom_tile() +
  scale_fill_manual(values=c("#FF0000", "#FFCCCC", "#E0E0E0", "#BABABA")) + theme_bw()
dev.off()

#########
consensus_tissue <- fread("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/HPA/rna_tissue_consensus.tsv.gz")
all <- consensus_tissue %>% filter(`Gene name` %in% c("ATP11A", "DPP9", "NPNT", "MUC1", "OAS1"))

library("corrplot")

all <- all %>% arrange(Tissue)
tmp <- all[,c("Gene name", "Tissue","Level", "nTPM")]
library(forcats)
# tmp <- tmp %>% mutate(level = fct_relevel(Level, 
#                                           "High", "Medium", "Low", 
#                                           "Not detected"),
#                       tissue = fct_relevel(Tissue, 
#                                            "lung", unique(all$Tissue)[-22]),
#                       genename = fct_relevel(`Gene name`, 
#                                              "ATP11A","DPP9", "NPNT", "OAS1", "MUC1"))
tmp <- tmp %>% mutate(tissue = fct_relevel(Tissue, 
                                                   c(unique(tmp$Tissue)[-24], "lung")))

tmp <- tmp %>% group_by(`Gene name`) %>% 
  mutate(`Normalized nTPM` = scale(nTPM))
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/Fig4A.png",width=400, height=800)
ggplot(tmp, aes(x=`Gene name`, y=tissue, fill=`Normalized nTPM`)) + geom_tile() + theme_bw() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab" )
dev.off()


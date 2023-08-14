setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/scRNAseq/Nature")

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
lung <- fread("GSE171668_lung_metadata.csv.gz")
lung <- lung %>% filter(doublet == FALSE)
celltypes <- unique(lung$manual_annotation_SubCluster)
library(Seurat)
genelist <- fread("../genelist_sQTL", header=F)
GSM5229973 <- Read10X_h5("GSM5229973_02-P005175-S053-R01_raw_feature_bc_matrix.h5")
GSM5230008 <- Read10X_h5("GSM5230008_04-P054921-S024-R01_raw_feature_bc_matrix.h5")
GSM5230014 <- Read10X_h5("GSM5230014_04-P079042-S015-R01_raw_feature_bc_matrix.h5")
GSM5230022 <- Read10X_h5("GSM5230022_04-P103142-S089-R01_raw_feature_bc_matrix.h5")
GSM5230023 <- Read10X_h5("GSM5230023_04-P103142-S109-R01_raw_feature_bc_matrix.h5")
GSM5230026 <- Read10X_h5("GSM5230026_04-P103142-S148-R01_raw_feature_bc_matrix.h5")
GSM5230027 <- Read10X_h5("GSM5230027_04-P103142-S149-R01_raw_feature_bc_matrix.h5")
GSM5230028 <- Read10X_h5("GSM5230028_04-P103142-S150-R01_raw_feature_bc_matrix.h5")
GSM5230030 <- Read10X_h5("GSM5230030_12-P230638-S005-R01_raw_feature_bc_matrix.h5")
GSM5230031 <- Read10X_h5("GSM5230031_12-P485759-S008-R01_raw_feature_bc_matrix.h5")
GSM5230037 <- Read10X_h5("GSM5230037_12-P617758-S008-R01_raw_feature_bc_matrix.h5")
GSM5230040 <- Read10X_h5("GSM5230040_12-P852049-S007-R01_raw_feature_bc_matrix.h5")
GSM5230044 <- Read10X_h5("GSM5230044_12-P890292-S007-R01_raw_feature_bc_matrix.h5")
GSM5230059 <- Read10X_h5("GSM5230059_02-P118946-S058-R01_out_FPR_0.01_filtered.h5")
GSM5230063 <- Read10X_h5("GSM5230063_02-P166169-S099-R01_out_FPR_0.01_filtered.h5")
GSM5230064 <- Read10X_h5("GSM5230064_02-P166169-S099-R02_out_FPR_0.01_filtered.h5")
GSM5230066 <- Read10X_h5("GSM5230066_02-P166169-S113-R01_out_FPR_0.01_filtered.h5")
GSM5230068 <- Read10X_h5("GSM5230068_02-P240970-S008-R01_out_FPR_0.01_filtered.h5")
GSM5230075 <- Read10X_h5("GSM5230075_02-P248880-S086-R01_out_FPR_0.01_filtered.h5")
GSM5230080 <- Read10X_h5("GSM5230080_02-P334354-S087-R01_out_FPR_0.01_filtered.h5")
GSM5230085 <- Read10X_h5("GSM5230085_02-P348762-S056-R01_out_FPR_0.01_filtered.h5")
GSM5230087 <- Read10X_h5("GSM5230087_04-P006354-S015-R01_out_FPR_0.01_filtered.h5")
GSM5230090 <- Read10X_h5("GSM5230090_04-P006354-S057-R01_out_FPR_0.01_filtered.h5")
GSM5230091 <- Read10X_h5("GSM5230091_04-P006354-S057-R02_out_FPR_0.01_filtered.h5")

for(i in seq(1, length(celltypes))){
  tmp <- lung %>% filter(manual_annotation_SubCluster == celltypes[i])
  tmp <- tmp %>% mutate(seq = paste0(str_split(barcodes, pattern="-", simplify = T)[,5],"-1"))
  tmp <- tmp %>% mutate(ID = paste0(str_split(barcodes, pattern="-", simplify = T)[,1],"-",str_split(barcodes, pattern="-", simplify = T)[,2],"-",
                                    str_split(barcodes, pattern="-", simplify = T)[,3],"-", str_split(barcodes, pattern="-", simplify = T)[,4]))
  out <- data.frame(matrix(0, length(genelist$V1), 26))
  colnames(out) <- c("celltype", "gene", unique(lung$donor))
  out$celltype <- celltypes[i]
  out$gene <- genelist$V1
  for(j in seq(6, length(genelist$V1))){
    GSM5229973_tmp <- GSM5229973[which(rownames(GSM5229973) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5229973_tmp <- GSM5229973_tmp[which(names(GSM5229973_tmp) %in% tmp$seq[tmp$ID == "02-P005175-S053-R01"])]
    out$D1[j] <- mean(GSM5229973_tmp, na.rm=T)
    GSM5230008_tmp <- GSM5230008[which(rownames(GSM5230008) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230008_tmp <- GSM5230008_tmp[which(names(GSM5230008_tmp) %in% tmp$seq[tmp$ID == "04-P054921-S024-R01"])]
    out$D10[j] <- mean(GSM5230008_tmp, na.rm=T)
    GSM5230014_tmp <- GSM5230014[which(rownames(GSM5230014) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230014_tmp <- GSM5230014_tmp[which(names(GSM5230014_tmp) %in% tmp$seq[tmp$ID == "04-P079042-S015-R01"])]
    out$D11[j] <- mean(GSM5230014_tmp, na.rm=T)
    GSM5230022_tmp <- GSM5230022[which(rownames(GSM5230022) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230022_tmp <- GSM5230022_tmp[which(names(GSM5230022_tmp) %in% tmp$seq[tmp$ID == "04-P103142-S089-R01"])]
    out$D12_1[j] <- mean(GSM5230022_tmp, na.rm=T)
    GSM5230023_tmp <- GSM5230023[which(rownames(GSM5230023) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230023_tmp <- GSM5230023_tmp[which(names(GSM5230023_tmp) %in% tmp$seq[tmp$ID == "04-P103142-S109-R01"])]
    out$D12_2[j] <- mean(GSM5230023_tmp, na.rm=T)
    GSM5230026_tmp <- GSM5230026[which(rownames(GSM5230026) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230026_tmp <- GSM5230026_tmp[which(names(GSM5230026_tmp) %in% tmp$seq[tmp$ID == "04-P103142-S148-R01"])]
    out$D12_3[j] <- mean(GSM5230026_tmp, na.rm=T)
    GSM5230027_tmp <- GSM5230027[which(rownames(GSM5230027) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230027_tmp <- GSM5230027_tmp[which(names(GSM5230027_tmp) %in% tmp$seq[tmp$ID == "04-P103142-S149-R01"])]
    out$D12_4[j] <- mean(GSM5230027_tmp, na.rm=T)
    GSM5230028_tmp <- GSM5230028[which(rownames(GSM5230028) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230028_tmp <- GSM5230028_tmp[which(names(GSM5230028_tmp) %in% tmp$seq[tmp$ID == "04-P103142-S150-R01"])]
    out$D12_5[j] <- mean(GSM5230028_tmp, na.rm=T)
    GSM5230030_tmp <- GSM5230030[which(rownames(GSM5230030) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230030_tmp <- GSM5230030_tmp[which(names(GSM5230030_tmp) %in% tmp$seq[tmp$ID == "12-P230638-S005-R01"])]
    out$D13[j] <- mean(GSM5230030_tmp, na.rm=T)
    GSM5230031_tmp <- GSM5230031[which(rownames(GSM5230031) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230031_tmp <- GSM5230031_tmp[which(names(GSM5230031_tmp) %in% tmp$seq[tmp$ID == "12-P485759-S008-R01"])]
    out$D14[j] <- mean(GSM5230031_tmp, na.rm=T)
    GSM5230037_tmp <- GSM5230037[which(rownames(GSM5230037) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230037_tmp <- GSM5230037_tmp[which(names(GSM5230037_tmp) %in% tmp$seq[tmp$ID == "12-P617758-S008-R01"])]
    out$D15[j] <- mean(GSM5230037_tmp, na.rm=T)
    GSM5230040_tmp <- GSM5230040[which(rownames(GSM5230040) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230040_tmp <- GSM5230040_tmp[which(names(GSM5230040_tmp) %in% tmp$seq[tmp$ID == "12-P852049-S007-R01"])]
    out$D16[j] <- mean(GSM5230040_tmp, na.rm=T)
    GSM5230044_tmp <- GSM5230044[which(rownames(GSM5230044) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230044_tmp <- GSM5230044_tmp[which(names(GSM5230044_tmp) %in% tmp$seq[tmp$ID == "12-P890292-S007-R01"])]
    out$D17[j] <- mean(GSM5230044_tmp, na.rm=T)
    GSM5230059_tmp <- GSM5230059[which(rownames(GSM5230059) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230059_tmp <- GSM5230059_tmp[which(names(GSM5230059_tmp) %in% tmp$seq[tmp$ID == "02-P118946-S058-R01"])]
    out$D2[j] <- mean(GSM5230059_tmp, na.rm=T)
    GSM5230063_tmp <- GSM5230063[which(rownames(GSM5230063) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230063_tmp <- GSM5230063_tmp[which(names(GSM5230063_tmp) %in% tmp$seq[tmp$ID == "02-P166169-S099-R01"])]
    out$D3_1[j] <- mean(GSM5230063_tmp, na.rm=T)
    GSM5230064_tmp <- GSM5230064[which(rownames(GSM5230064) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230064_tmp <- GSM5230064_tmp[which(names(GSM5230064_tmp) %in% tmp$seq[tmp$ID == "02-P166169-S099-R02"])]
    out$D3_2[j] <- mean(GSM5230064_tmp, na.rm=T)
    GSM5230066_tmp <- GSM5230066[which(rownames(GSM5230066) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230066_tmp <- GSM5230066_tmp[which(names(GSM5230066_tmp) %in% tmp$seq[tmp$ID == "02-P166169-S113-R01"])]
    out$D3_2[j] <- mean(GSM5230066_tmp, na.rm=T)
    GSM5230068_tmp <- GSM5230068[which(rownames(GSM5230068) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230068_tmp <- GSM5230068_tmp[which(names(GSM5230068_tmp) %in% tmp$seq[tmp$ID == "02-P240970-S008-R01"])]
    out$D4[j] <- mean(GSM5230068_tmp, na.rm=T)
    GSM5230075_tmp <- GSM5230075[which(rownames(GSM5230075) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230075_tmp <- GSM5230075_tmp[which(names(GSM5230075_tmp) %in% tmp$seq[tmp$ID == "02-P248880-S086-R01"])]
    out$D5[j] <- mean(GSM5230075_tmp, na.rm=T)
    GSM5230080_tmp <- GSM5230080[which(rownames(GSM5230080) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230080_tmp <- GSM5230080_tmp[which(names(GSM5230080_tmp) %in% tmp$seq[tmp$ID == "02-P334354-S087-R01"])]
    out$D6[j] <- mean(GSM5230080_tmp, na.rm=T)
    GSM5230085_tmp <- GSM5230085[which(rownames(GSM5230085) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230085_tmp <- GSM5230085_tmp[which(names(GSM5230085_tmp) %in% tmp$seq[tmp$ID == "02-P348762-S056-R01"])]
    out$D7[j] <- mean(GSM5230085_tmp, na.rm=T)
    GSM5230087_tmp <- GSM5230087[which(rownames(GSM5230087) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230087_tmp <- GSM5230087_tmp[which(names(GSM5230087_tmp) %in% tmp$seq[tmp$ID == "04-P006354-S015-R01"])]
    out$D8_1[j] <- mean(GSM5230087_tmp, na.rm=T)
    GSM5230090_tmp <- GSM5230090[which(rownames(GSM5230090) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230090_tmp <- GSM5230090_tmp[which(names(GSM5230090_tmp) %in% tmp$seq[tmp$ID == "04-P006354-S057-R01"])]
    out$D8_2[j] <- mean(GSM5230090_tmp, na.rm=T)
    GSM5230091_tmp <- GSM5230091[which(rownames(GSM5230091) == paste0("GRCh38premrna_",genelist$V1[j])),]
    GSM5230091_tmp <- GSM5230091_tmp[which(names(GSM5230091_tmp) %in% tmp$seq[tmp$ID == "04-P006354-S057-R02"])]
    out$D8_3[j] <- mean(GSM5230091_tmp, na.rm=T)
  }
  out %>% write.table(file="mean_expression_per_celltype.tsv", quote=F, col.names = F, row.names = F, sep="\t", append=T)
}

out <- fread("mean_expression_per_celltype.tsv")
colnames(out) <- c("celltype", "gene", unique(lung$donor))
out %>% write.table(file="mean_expression_per_celltype.tsv", quote=F, col.names = T, row.names = F, sep="\t", append=F)

out$mean <- rowMeans(out[,c(3:26)], na.rm=T)

library(ggplot2)
library(ggridges)
theme_set(theme_minimal())

data <- melt(setDT(out), id.vars = c("celltype","gene"), 
             variable.name = "sample")

data <- data %>% filter(!grepl("Doublet|doublet", celltype))

png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/Fig4C.png",width=1000, height=400)
data %>% ggplot(aes(x = celltype, y = value, color=celltype)) +
  geom_jitter(size=0.1) + facet_wrap(~gene, ncol=1, scales = "free_y") + geom_boxplot(alpha=0.6,trim = FALSE) +
  theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        strip.text = element_text(size=10)
        )
dev.off()



library(tximport)
library(tidyverse)
library(data.table)
library(openxlsx)

batch <- read.xlsx("/project/richards/restricted/bqc19_release6/RNA_data/BQC19_RNAsample list_GQ.xlsx", sheet = "Sheet0")
# batch <- batch %>% mutate(SampleID = case_when(grepl("VAR", Nom) ~ Nom,
#                                                Extraction == "CHUM" ~ paste0("VAR", Nom)))
# 
batch <- batch %>% mutate(Run = `Run.#`,
                          SampleID = Nom)
batch <- batch %>% dplyr::select(SampleID, Run, Extraction, Date.QC, Date.Run)
batch <- batch %>% mutate(Date.QC = as.Date(Date.QC, format="%d/%m/%Y"),
                          Date.Run = as.Date(Date.Run, origin="1899-12-30"))
batch <- batch %>% arrange(desc(Date.Run))
batch <- batch[!duplicated(batch$SampleID),]

classfication <- fread("/scratch/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/phenotype/release6/bqc_blood_samples_with_Inf_nonInf_classification_with_severity.tsv")
classfication <- classfication %>% dplyr::select(BQCID, EventName, RepeatInstance, covid19_test_result, first_covid19_test_date, sx_date, sample_status, priority,
                                                 infectiousWithLongTerm, A2, B2, C2, final_severity)
mapping <- fread("/scratch/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/phenotype/release6/paxsample.mapping.tsv")
mapping <- mapping %>% inner_join(classfication, by=c("BQCID"="BQCID", "redcap_event_name"="EventName",
                                                      "redcap_repeat_instance"="RepeatInstance"))

#remove duplicated vcodes
mapping <- mapping %>% filter(!(vcode %in% mapping$vcode[duplicated(mapping$vcode)]))
batch <- batch %>% inner_join(mapping, by=c("SampleID"="vcode"))

kallisto.gene <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/05.BQC/RNAseq_Kallisto/Kallisto.gene.Precombat.rds")
kallisto.gene.postCombat <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/05.BQC/RNAseq_Kallisto/Kallisto.gene.combat.rds")

library(edgeR)
edata_preCombat <- as.matrix(kallisto.gene)
edata_postCombat <- as.matrix(kallisto.gene.postCombat)

cor.test(edata_preCombat, edata_postCombat)

edata_preCombat.pca <- prcomp(t(edata_preCombat), scale = TRUE)
edata_postCombat.pca <- prcomp(t(edata_postCombat), scale = TRUE)

predat.umap <- uwot::umap(edata_preCombat.pca$x, verbose = T)
postdat.umap <- uwot::umap(edata_postCombat.pca$x, verbose = T)

sample3 <- readRDS("Kallisto.batch.rds")

tmp <- data.frame(edata_preCombat.pca$x)
tmp1 <- data.frame(predat.umap[,1:2])
colnames(tmp1) <- c("UMAP1", "UMAP2")

final <- bind_cols(sample3, tmp, tmp1)

ggplot(final, aes(x=PC1, y=PC2, color=as.factor(final_severity))) + geom_point()
ggplot(final, aes(x=PC3, y=PC4, color=as.factor(final_severity))) + geom_point()
ggplot(final, aes(x=PC5, y=PC6, color=as.factor(final_severity))) + geom_point()
ggplot(final, aes(x=PC7, y=PC8, color=as.factor(final_severity))) + geom_point()

ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(final_severity))) + geom_point()


ggplot(final, aes(x=PC1, y=PC2, color=as.factor(samplestatus))) + geom_point()
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(samplestatus))) + geom_point()

ggplot(final, aes(x=PC1, y=PC2, color=as.factor(batch))) + geom_point()
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(Extraction))) + geom_point()


sample3 <- readRDS("Kallisto.batch.rds")

tmp <- data.frame(edata_postCombat.pca$x)
tmp1 <- data.frame(postdat.umap[,1:2])
colnames(tmp1) <- c("UMAP1", "UMAP2")

final <- bind_cols(sample3, tmp, tmp1)

ggplot(final, aes(x=PC1, y=PC2, color=as.factor(final_severity))) + geom_point()
ggplot(final, aes(x=PC3, y=PC4, color=as.factor(final_severity))) + geom_point()
ggplot(final, aes(x=PC5, y=PC6, color=as.factor(final_severity))) + geom_point()
ggplot(final, aes(x=PC7, y=PC8, color=as.factor(final_severity))) + geom_point()

ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(final_severity))) + geom_point()


ggplot(final, aes(x=PC1, y=PC2, color=as.factor(samplestatus))) + geom_point()
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(samplestatus))) + geom_point()

ggplot(final, aes(x=PC3, y=PC4, color=as.factor(Extraction))) + geom_point()
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(Run))) + geom_point()

##GTEx

GTEx_cov <- fread("/project/richards/restricted/bqc19_release6/RNA_data/results/data_mirror_sample_set/release6.3_2022-05-05_inf_all_pop_all/expression_combined_covariates/release6.3_2022-05-05_inf_all_pop_all.combined_covariates.txt")
GTEx_data_ori <- fread("/project/richards/restricted/bqc19_release6/RNA_data/results/data_mirror_sample_set/release6.3_2022-05-05_inf_all_pop_all/expression_bed/release6.3_2022-05-05_inf_all_pop_all.expression.bed.gz")

GTEx_data <- GTEx_data_ori[,c(-1,-2,-3,-4)]
GTEx_data <- data.frame(t(GTEx_data))
colnames(GTEx_data) <- GTEx_data_ori$gene_id
GTEx_data$BQCID <- colnames(GTEx_data_ori)[c(-1,-2,-3,-4)]

library("jsonlite")
mapping <- fromJSON("/project/richards/restricted/bqc19_release6/RNA_data/results/mirror.sample_set.json", simplifyVector = TRUE)
all_sample <- mapping$attributes$sample[,2][[2]]
all_sample <- all_sample %>% mutate(BQCID = str_split(entityName, pattern="\\.", simplify = T)[,1],
                                    SampleID = str_split(entityName, pattern="\\.", simplify = T)[,2])
all_sample <- all_sample %>% dplyr::select(BQCID, SampleID)

GTEx_data <- GTEx_data %>% inner_join(all_sample, by="BQCID")

GTEx_data <- GTEx_data %>% mutate(pipeline="gtex")
colnames(GTEx_data) <- str_split(colnames(GTEx_data), pattern="\\.", simplify = T)[,1]
Kallisto_comp <- data.frame(t(as.matrix(kallisto.gene)))

Kallisto_comp <- Kallisto_comp %>% mutate(pipeline = "kallisto")
colnames(Kallisto_comp) <- str_split(colnames(Kallisto_comp), pattern="\\.", simplify = T)[,1]
Kallisto <- bind_cols(sample3, Kallisto_comp)

final <- bind_rows(GTEx_data, Kallisto) %>% drop_na(SampleID)
final <- final %>% group_by(SampleID) %>%
  filter(length(SampleID) == 2)

final <- final %>% 
  dplyr::select(
    where(
      ~!any(is.na(.x))
    )
  )

kallisto <- final %>% filter(pipeline == "kallisto") %>% arrange(SampleID) %>% dplyr::select(colnames(final)[grepl("ENSG", colnames(final))])
gtex <- final %>% filter(pipeline == "gtex") %>% arrange(SampleID) %>% dplyr::select(colnames(final)[grepl("ENSG", colnames(final))])

cor.test(as.matrix(kallisto[,-1]), as.matrix(gtex[,-1]))
#0.2818946


###GTEx expression
GTEx_data <- GTEx_data_ori[,c(-1,-2,-3,-4)]
GTEx_data <- data.frame(t(GTEx_data))
colnames(GTEx_data) <- GTEx_data_ori$gene_id
GTEx_data$BQCID <- colnames(GTEx_data_ori)[c(-1,-2,-3,-4)]

library("jsonlite")
mapping <- fromJSON("/project/richards/restricted/bqc19_release6/RNA_data/results/mirror.sample_set.json", simplifyVector = TRUE)
all_sample <- mapping$attributes$sample[,2][[2]]
all_sample <- all_sample %>% mutate(BQCID = str_split(entityName, pattern="\\.", simplify = T)[,1],
                                    SampleID = str_split(entityName, pattern="\\.", simplify = T)[,2])
all_sample <- all_sample %>% dplyr::select(BQCID, SampleID)

GTEx_data <- GTEx_data %>% inner_join(all_sample, by="BQCID")

batch <- read.xlsx("/project/richards/restricted/bqc19_release6/RNA_data/BQC19_RNAsample list_GQ.xlsx", sheet = "Sheet0")
# batch <- batch %>% mutate(SampleID = case_when(grepl("VAR", Nom) ~ Nom,
#                                                Extraction == "CHUM" ~ paste0("VAR", Nom)))
# 
batch <- batch %>% mutate(Run = `Run.#`,
                          SampleID = Nom)
batch <- batch %>% dplyr::select(SampleID, Run, Extraction, Date.QC, Date.Run)
batch <- batch %>% mutate(Date.QC = as.Date(Date.QC, format="%d/%m/%Y"),
                          Date.Run = as.Date(Date.Run, origin="1899-12-30"))
batch <- batch %>% arrange(desc(Date.Run))
batch <- batch[!duplicated(batch$SampleID),]

#sample <- sample %>% dplyr::rename(sampleID = `entity:sample_id`)

batch <- batch %>%  dplyr::filter(SampleID %in% sample$VAR)

classfication <- fread("/scratch/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/phenotype/release6/bqc_blood_samples_with_Inf_nonInf_classification_with_severity.tsv")
classfication <- classfication %>% dplyr::select(BQCID, EventName, RepeatInstance, covid19_test_result, first_covid19_test_date, sx_date, sample_status, priority,
                                                 infectiousWithLongTerm, A2, B2, C2, final_severity)
mapping <- fread("/scratch/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/phenotype/release6/paxsample.mapping.tsv")
mapping <- mapping %>% inner_join(classfication, by=c("BQCID"="BQCID", "redcap_event_name"="EventName",
                                                      "redcap_repeat_instance"="RepeatInstance"))

#remove duplicated vcodes
mapping <- mapping %>% dplyr::filter(!(vcode %in% mapping$vcode[duplicated(mapping$vcode)]))


batch <- batch %>% inner_join(mapping, by=c("SampleID"="vcode"))

GTEx_data <- GTEx_data %>% inner_join(batch, by=c("BQCID"="BQCID", "SampleID"="SampleID"))

GTEx_data.pca <- prcomp(GTEx_data[,colnames(GTEx_data)[grepl("ENSG", colnames(GTEx_data))]], scale = TRUE)

GTEx_data.umap <- uwot::umap(GTEx_data.pca$x, verbose = T)

tmp <- data.frame(GTEx_data.pca$x)
tmp1 <- data.frame(GTEx_data.umap[,1:2])
colnames(tmp1) <- c("UMAP1", "UMAP2")

final <- bind_cols(GTEx_data[,colnames(GTEx_data)[!grepl("ENSG", colnames(GTEx_data))]], tmp, tmp1)

library(cowplot)
p1 <- ggplot(final, aes(x=PC1, y=PC2, color=as.factor(sample_status))) + geom_point()
p2 <- ggplot(final, aes(x=PC3, y=PC4, color=as.factor(sample_status))) + geom_point()
p3 <- ggplot(final, aes(x=PC5, y=PC6, color=as.factor(sample_status))) + geom_point()
p4 <- ggplot(final, aes(x=PC7, y=PC8, color=as.factor(sample_status))) + geom_point()
plot_grid(p1, p2, p3, p4)
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(final_severity))) + geom_point()
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(Run))) + geom_point()
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(Extraction))) + geom_point()
ggplot(final, aes(x=UMAP1, y=UMAP2, color=as.factor(sample_status))) + geom_point()

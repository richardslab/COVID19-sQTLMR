setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/05.BQC/RNAseq_Kallisto")

library(ggplot2)
library(tximport)
library(tidyverse)
library(data.table)
library(openxlsx)
library(stringr)

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


##transcript
kallisto.transcript.postCombat <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/05.BQC/RNAseq_Kallisto/Kallisto.transcript.combat.rds")
data <- data.frame(t(kallisto.transcript.postCombat))
data$SampleID <- rownames(data)
colnames(data) <- str_split(colnames(data), pattern = "\\.", simplify = T)[,1]

##transcriptlist 
genelist <- c("ENSG00000168743", "ENSG00000175164", "ENSG00000068650", "ENSG00000160783", "ENSG00000089127",
              "ENSG00000160766", "ENSG00000169231", "ENSG00000142002", "ENSG00000185303", "ENSG00000185499", "ENSG00000091592")
require(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "gene_biotype"),
  uniqueRows=TRUE)
annotLookup <- annotLookup %>% dplyr::filter(ensembl_gene_id %in% genelist)

data <- data %>% dplyr::select(SampleID, any_of(annotLookup$ensembl_transcript_id))

sample3 <- readRDS("Kallisto.batch.rds")

pheno_raw <- read.csv("/project/richards/restricted/bqc19_release6/Clinical_Phenotypic/redcap_clinical_data_raw_2022-03-24.csv", header = T, fileEncoding="latin1", sep="\t")
pheno <- pheno_raw %>% dplyr::select(BQCID, age, female) %>% drop_na()

data <- data %>% inner_join(sample3, by="SampleID")
data <- data %>% left_join(pheno, by="BQCID")

data <- data %>% mutate(time = as.Date(draw_date) - as.Date(sx_date))
##trajectory
png("tmp.png", width=700, height = 400)
ggplot(data, aes(x=time, y=ENST00000262960, col=final_severity)) + geom_line(alpha=0.2,size=0.4,aes(group=BQCID,col=final_severity))+ 
  theme_bw() + geom_smooth(se=T,size=1, aes(group=final_severity,col=final_severity)) +
  xlab("Days from symptoms") + 
  ylab("") +
  theme(legend.title = element_text( size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  scale_colour_brewer(palette="Dark2") +
  guides(colour = guide_legend(reverse=TRUE))
dev.off()

##infectious
inf <- data %>% group_by(BQCID) %>%
  dplyr::filter(infectiousWithLongTerm == "Yes") %>%
  dplyr::filter(priority == max(priority, na.rm=T))
transcriptlist <- colnames(data)[2:21]
  
out <- data.frame(matrix(0, length(transcriptlist), 4))
colnames(out) <- c("transcript", "beta", "se", "pval")

for(i in seq(1,length(transcriptlist))){
  gene <- transcriptlist[i]
  LM <- glm(paste0("A2 ~",gene," + age + female"), data=inf, family="binomial")
  out$transcript[i] <- gene
  out[i,2:4] <- summary(LM)$coefficients[2,c(1,2,4)]
}

tmp <- annotLookup %>% filter(hgnc_symbol == "DPP9") %>% select(ensembl_transcript_id)
inf %>% select(any_of(tmp$ensembl_transcript_id))

inf <- inf %>% mutate(ENSG00000142002 = ENST00000262960 + ENST00000597145)

LM <- glm(paste0("A2 ~ ENST00000597145 + age + female"), data=inf, family="binomial")
##noninfectious
noninf <- data %>% group_by(BQCID) %>%
  dplyr::filter(infectiousWithLongTerm == "No") %>%
  dplyr::filter(priority == max(priority, na.rm=T))

out <- data.frame(matrix(0, length(transcriptlist), 4))
colnames(out) <- c("transcript", "beta", "se", "pval")

for(i in seq(1,length(transcriptlist))){
  gene <- transcriptlist[i]
  LM <- glm(paste0("A2 ~",gene," + age + female"), data=noninf, family="binomial")
  out$transcript[i] <- gene
  out[i,2:4] <- summary(LM)$coefficients[2,c(1,2,4)]
}

noninf <- noninf %>% dplyr::filter(final_severity == "0NonCOVID")

tmp <- bind_rows(inf, noninf)

LM <- glm(paste0("A2 ~ ENST00000375645 + age + female"), data=tmp, family="binomial")
summary(LM)
###gene
kallisto.gene.postCombat <- readRDS("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/05.BQC/RNAseq_Kallisto/Kallisto.gene.combat.rds")
data <- data.frame(t(kallisto.gene.postCombat))
data$SampleID <- rownames(data)

##genelist 
genelist <- c("ENSG00000168743", "ENSG00000175164", "ENSG00000068650", "ENSG00000160783", "ENSG00000089127",
              "ENSG00000160766", "ENSG00000169231", "ENSG00000142002", "ENSG00000185303", "ENSG00000185499", "ENSG00000091592",
              "ENSG00000105483")
data <- data %>% dplyr::select(SampleID, any_of(genelist))

sample3 <- readRDS("Kallisto.batch.rds")

pheno_raw <- read.csv("/project/richards/restricted/bqc19_release6/Clinical_Phenotypic/redcap_clinical_data_raw_2022-03-24.csv", header = T, fileEncoding="latin1", sep="\t")
pheno <- pheno_raw %>% dplyr::select(BQCID, age, female) %>% drop_na()

data <- data %>% inner_join(sample3, by="SampleID")
data <- data %>% left_join(pheno, by="BQCID")

data <- data %>% mutate(time = as.Date(draw_date) - as.Date(sx_date))
##trajectory
png("tmp.png", width=700, height = 400)
ggplot(data, aes(x=time, y=ENSG00000089127, col=final_severity)) + geom_line(alpha=0.2,size=0.4,aes(group=BQCID,col=final_severity))+ 
  theme_bw() + geom_smooth(se=T,size=1, aes(group=final_severity,col=final_severity)) +
  xlab("Days from symptoms") + 
  ylab("") +
  theme(legend.title = element_text( size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 20)) +
  scale_colour_brewer(palette="Dark2") +
  guides(colour = guide_legend(reverse=TRUE))
dev.off()

##infectious
inf <- data %>% group_by(BQCID) %>%
  dplyr::filter(infectiousWithLongTerm == "Yes") %>%
  dplyr:: filter(priority == max(priority, na.rm=T))
genelist <- c("ENSG00000068650", "ENSG00000160783", "ENSG00000089127", "ENSG00000169231", "ENSG00000142002","ENSG00000091592","ENSG00000105483")

out <- data.frame(matrix(0, length(genelist), 4))
colnames(out) <- c("gene", "beta", "se", "pval")

for(i in seq(1,length(genelist))){
  gene <- genelist[i]
  LM <- glm(paste0("A2 ~",gene," + age + female"), data=inf, family="binomial")
  out$gene[i] <- gene
  out[i,2:4] <- summary(LM)$coefficients[2,c(1,2,4)]
}

LM <- glm(paste0("A2 ~ ENSG00000142002 + age + female"), data=inf, family="binomial")
summary(LM)


##noninfectious
noninf <- data %>% group_by(BQCID) %>%
  filter(infectiousWithLongTerm == "No") %>%
  filter(priority == max(priority, na.rm=T))

out <- data.frame(matrix(0, length(genelist), 4))
colnames(out) <- c("gene", "beta", "se", "pval")
for(i in seq(1,length(genelist))){
  gene <- genelist[i]
  LM <- glm(paste0("A2 ~",gene," + age + female"), data=noninf, family="binomial")
  out$gene[i] <- gene
  out[i,2:4] <- summary(LM)$coefficients[2,c(1,2,4)]
}

LM <- glm(paste0("A2 ~ ENSG00000142002 + age + female"), data=noninf, family="binomial")
summary(LM)



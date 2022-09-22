setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/data/05.BQC/RNAseq_Kallisto")
library(tximport)
library(tidyverse)
library(data.table)
library(openxlsx)
kallisto.files <- file.path(list.files('.', pattern = 'BQC'), 'abundance.h5')
kallisto.files <- kallisto.files[-1079]
Names <- kallisto.files
Names <- data.frame(Names)
Names <- Names %>% mutate(Names = stringr::str_split(Names, pattern = "/", simplify = T)[,1])
Names <- Names %>% mutate(Names = stringr::str_split(Names, pattern = "\\.", simplify = T)[,2])

names(kallisto.files) <- Names$Names
tx.exp <- tximport(kallisto.files, type = "kallisto", txOut = TRUE)
saveRDS(tx.exp, file="Kallisto.transcript.rds")
tx.exp <- readRDS(file="Kallisto.transcript.rds")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
k <- keys(edb, keytype = "TXID")
tx2gene <- select(edb, k, "GENEID", "TXID")

tx <- data.frame(row.names(tx.exp$abundance))
names(tx) <- c("TXID")
tx <- tx %>% mutate(key = str_split(TXID, pattern ="\\.", simplify = T)[,1])
tx <- tx %>% left_join(tx2gene, by=c("key"="TXID"), sort=F)

tx2gene <- tx %>% dplyr::select(TXID, GENEID)
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene <- tx2gene %>% drop_na(GENEID)
saveRDS(tx2gene, file="tx2gene.rds")

tmp <- kallisto.files
txi.kallisto.tsv <- tximport(tmp, type = "kallisto", tx2gene = tx2gene,
                             txOut = FALSE)
data <- txi.kallisto.tsv$counts
d0 <- DGEList(data)

saveRDS(txi.kallisto.tsv, file="Kallisto.gene.rds")


tx.exp <- readRDS(file="Kallisto.transcript.rds")
# BiocManager::install("edgeR")
# helpful sourse https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

library(edgeR)
data <- tx.exp$counts
d0 <- DGEList(data)

##sample QC done by Yossi
sample <- fread("sample.tsv")
sample <- sample %>% dplyr::filter(is.na(blacklisted_Q) | blacklisted_Q != 1)
sample <- sample %>% mutate(VAR = str_split(`entity:sample_id`,pattern="\\.", simplify = T)[,2])
#batch
# library(readxlsb)
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

d1 <- d0[,c(batch$SampleID)]

# RNA filter1 median logCPM > 1
d1 <- calcNormFactors(d1)
cutoff <- 1
drop <- which(apply(log2(cpm(d1)), 1, median) <= cutoff)
dim(d1)
d2 <- d1[-drop,] 
dim(d2) # number of genes left

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
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

annotLookup <- annotLookup %>% dplyr::filter(gene_biotype == "protein_coding")

keep <- data.frame(rownames(d2))
keep <- keep %>% mutate(transcript = str_split(rownames.d2., pattern = "\\.", simplify = T)[,1])
keep <- keep %>% dplyr::filter(transcript %in% annotLookup$ensembl_transcript_id)
d2 <- d2[keep$rownames.d2.,]
dim(d2)
#convert to logCPM value
d2 <- calcNormFactors(d2)
d3 <- voom(d2)

saveRDS(d3, "Kallisto.transcript.Precombat.rds")
library(sva)
sample3 <- data.frame(colnames(d3))
colnames(sample3) <- "SampleID"
sample3 <- inner_join(sample3, batch, by="SampleID", sort = F)

sample3 <- sample3 %>% mutate(samplestatus = case_when(grepl("^Infectious", sample_status) ~ "Infectious",
                                                       grepl("^NonInfectious", sample_status) ~ "NonInfectious",
                                                       grepl("^LongTerm", sample_status) ~ "NonInfectious",
                                                       TRUE ~ "Others"))
sample3 <- sample3 %>% mutate(group=paste0(samplestatus,":",final_severity))
sample3 <- sample3 %>% mutate(batch = paste0(Run,":",Extraction))
edata <- as.matrix(d3)

mod = model.matrix(~as.factor(group), data=sample3)

group <- sample3$Extraction

d4 <- ComBat(as.matrix(d3), group, mod, par.prior=TRUE, prior.plots=FALSE)

saveRDS(d4, "Kallisto.transcript.combat.rds")


tx.exp <- readRDS(file="Kallisto.gene.rds")
# BiocManager::install("edgeR")
# helpful sourse https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

library(edgeR)
data <- tx.exp$counts
d0 <- DGEList(data)

d1 <- d0[,c(batch$SampleID)]

# RNA filter1 median logCPM > 1
d1 <- calcNormFactors(d1)
cutoff <- 1
drop <- which(apply(log2(cpm(d1)), 1, median) <= cutoff)
dim(d1)
d2 <- d1[-drop,] 
dim(d2) # number of genes left

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
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

annotLookup <- annotLookup %>% dplyr::filter(gene_biotype == "protein_coding")

keep <- data.frame(rownames(d2))
colnames(keep) <- "gene"
keep <- keep %>% dplyr::filter(gene %in% annotLookup$ensembl_gene_id)
d2 <- d2[keep$gene,]
dim(d2)
#convert to logCPM value
d2 <- calcNormFactors(d2)
d3 <- voom(d2)
saveRDS(d3, "Kallisto.gene.Precombat.rds")

library(sva)
sample3 <- data.frame(colnames(d3))
colnames(sample3) <- "SampleID"
sample3 <- inner_join(sample3, batch, by="SampleID", sort = F)

sample3 <- sample3 %>% mutate(samplestatus = case_when(grepl("^Infectious", sample_status) ~ "Infectious",
                                                       grepl("^NonInfectious", sample_status) ~ "NonInfectious",
                                                       grepl("^LongTerm", sample_status) ~ "NonInfectious",
                                                       TRUE ~ "Others"))
sample3 <- sample3 %>% mutate(group=paste0(samplestatus,":",final_severity))
sample3 <- sample3 %>% mutate(batch = paste0(Run,":",Extraction))

sample3 %>% saveRDS("Kallisto.batch.rds")
edata <- as.matrix(d3)

mod = model.matrix(~as.factor(group), data=sample3)

#Extraction correct batch effect than using batch..
group <- sample3$Extraction

d4 <- ComBat(as.matrix(d3), group, mod, par.prior=TRUE, prior.plots=FALSE)

saveRDS(d4, "Kallisto.gene.combat.rds")

rownames(d4)[which(rownames(d4) == "ENSG00000169231")]


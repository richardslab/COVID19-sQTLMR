
map1 <- fread("/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/sampleid.studyid.sex.map", header=F)
map1 <- map1 %>% drop_na(V1)
map2 <- read.xlsx("/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/Master_list_WHOLE_BLOOD_CORRECTED_10112021.xlsx", sheet = "concilie")
map2 <- map2 %>% dplyr::select(Alias.BQCid, Collection.site, Individual.id)
map2 <- map2 %>% drop_na(Alias.BQCid)
map2 <- map2 %>% unique()
map2 <- map2  %>% left_join(map1, by=c("Individual.id"="V2"))

colnames(map2) <- c("Alias.BQCid", "Collection.site", "Individual.id","genotypeID", "SEX")

write.table(map2, file="/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/sampleid.studyid.sex.map.BQCID",
            quote=F, col.names = T, row.names = F, sep="\t")

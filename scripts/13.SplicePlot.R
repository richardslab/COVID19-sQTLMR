setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/14.COVID19-esQTLMR/")

library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(Gviz)

edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "12")

gnm <- GRanges("12:112917700-112919389")

OAS1 <- proteins(edbx, filter = ~ genename == "OAS1")
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx)
gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000452357","ENST00000202917"),]
gnm_gns$symbol <- gnm_gns$transcript
gtx <- GeneRegionTrack(gnm_gns, geneSymbol = FALSE,
                       showId = TRUE)

## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm, fontsize=12)
## plot the region
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/ExtendedDataFig2A.png",width=800, height=400)
plotTracks(list(ht))
dev.off()

edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "13")

gnm_tx <- genomeToTranscript(gnm, edbx)

gnm <- GRanges("13:112875941-112881427")

ATP11A <- proteins(edbx, filter = ~ genename == "ATP11A")
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx, filter = GRangesFilter(gnm))
gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000375645","ENST00000375630"),]
gnm_gns$symbol <- gnm_gns$transcript
gtx <- GeneRegionTrack(gnm_gns, geneSymbol = FALSE,
                       showId = TRUE)

## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm, fontsize=12)
## plot the region
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/ExtendedDataFig2B.png",width=800, height=400)
plotTracks(list(ht))
dev.off()

#grepl(substr(ATP11A$protein_sequence[7], 1110, 1134), ATP11A$protein_sequence[8])
#1110-

edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "19")
gnm <- GRanges("19:4714337-4717615")

DPP9 <- proteins(edbx, filter = ~ genename == "DPP9")
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx)
# gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000262960","ENST00000601130","ENST00000600621","ENST00000594671","ENST00000598800","ENST00000703254",
# "ENST00000597849","ENST00000646573","ENST00000600556","ENST00000599163","ENST00000601720","ENST00000597900","ENST00000601173","ENST00000598360","ENST00000600497",
# "ENST00000597726","ENST00000602161","ENST00000593973","ENST00000595940","ENST00000598041","ENST00000597024","ENST00000595327","ENST00000599248",
# "ENST00000597145","ENST00000703255","ENST00000601764","ENST00000703253","ENST00000597253","ENST00000599998"),]

gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000598800", "ENST00000262960"),]

gnm_gns$symbol <- gnm_gns$transcript
gtx <- GeneRegionTrack(gnm_gns, geneSymbol = FALSE,
                       showId = TRUE)


## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm, fontsize=20)
## plot the region
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/ExtendedDataFig2C.png",width=800, height=300)

plotTracks(list(ht))
dev.off()
 
#ENST00000375645 nchar(ATP11A$protein_sequence[7])#1134
#ENST00000375630 nchar(ATP11A$protein_sequence[8])#1191

#ENST00000262960 nchar(DPP9$protein_sequence[2])
#ENST00000598800 nchar(DPP9$protein_sequence[4])

grepl(substr(DPP9$protein_sequence[4], 1, 863), DPP9$protein_sequence[2])


edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "4")
gnm <- GRanges("4:105898001-105927336")

NPNT <- proteins(edbx, filter = ~ genename == "NPNT")
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx)
gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000379987","ENST00000504787","ENST00000503451"),]

# gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000379987","ENST00000305572","ENST00000427316","ENST00000453617",
# "ENST00000514622","ENST00000506666","ENST00000503451","ENST00000514837","ENST00000504787","ENST00000504304","ENST00000505917",
#                                              "ENST00000513430","ENST00000505821","ENST00000506056","ENST00000514632","ENST00000511518"),]
gnm_gns$symbol <- gnm_gns$transcript
gtx <- GeneRegionTrack(gnm_gns, geneSymbol = FALSE,
                       showId = TRUE)

## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm, fontsize=20)
## plot the region
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/ExtendedDataFig2D.png",width=800, height=300)

plotTracks(list(ht))
dev.off()


edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "1")
gnm <- GRanges("1:155192310-155192786")

MUC1 <- proteins(edbx, filter = ~ genename == "MUC1")
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx)
gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000620103","ENST00000485118","ENST00000615517"),]
# gnm_gns <- gnm_gns[gnm_gns$transcript %in% c("ENST00000620103","ENST00000457295","ENST00000462317","ENST00000338684",
#                                              "ENST00000610359","ENST00000368392","ENST00000615517","ENST00000368393",
#                                              "ENST00000438413","ENST00000337604","ENST00000368390","ENST00000368398",
#                                              "ENST00000485118","ENST00000343256","ENST00000471283","ENST00000368389",
#                                              "ENST00000342482","ENST00000368396"),]
gnm_gns$symbol <- gnm_gns$transcript
gtx <- GeneRegionTrack(gnm_gns, geneSymbol = FALSE,
                       showId = TRUE)


## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm, fontsize=20)
## plot the region
png("/home/richards/tomoko.nakanishi/my_project/repo/COVID19-sQTLMR/results/ExtendedDataFig2E.png",width=800, height=500)
plotTracks(list(ht))
dev.off()

#ENST00000375645 nchar(ATP11A$protein_sequence[7])#1134
#ENST00000375630 nchar(ATP11A$protein_sequence[8])#1191

#ENST00000262960 nchar(DPP9$protein_sequence[2])
#ENST00000598800 nchar(DPP9$protein_sequence[4])

grepl(substr(DPP9$protein_sequence[4], 1, 863), DPP9$protein_sequence[2])




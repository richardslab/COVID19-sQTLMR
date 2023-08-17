library(tidyverse) ; library(vroom) ; library(EnvStats) ; library(vcfR) ; 
library(glue) ; library(data.table) ; library(magrittr) ; library(stringr)
library(rtracklayer)

setwd("/scratch/richards/julian.willett/14.GTEx_sQTL_Tomoko_Work")
source('Functions.R')

# imput vcf file has skip=15

vars = c('chr19:4714337:4717615:clu_25374:ENSG00000142002.16', #GRCh38
         'chr13:112875941:112880546:clu_3196:ENSG00000068650.18',
         'chr1:155192310:155192786:clu_59160:ENSG00000185499.16',
         'chr12:112917700:112919389:clu_6211:ENSG00000089127.12',
         'chr4:105898001:105927336:clu_44968:ENSG00000168743.12',
         'chr1:156233728:156236349:clu_59301:ENSG00000160783.19')
plot.titles = c('MUC1 chr1:155192310:155192786',
                'PMF1 chr1:156233728:156236349',
                'NPNT chr4:105898001:105927336',
                'OAS1 chr12:112917700:112919389',
                'ATP11A chr13:112875941:112880546',
                'DPP9 chr19:4714337:4717615') #all grch38

quant.excised = vroom('Lung.v8.leafcutter_phenotypes.bed.gz') %>% #GRCh38
  dplyr::select(-start,-end) %>%
  dplyr::filter(`#Chr` %in% c('chr1','chr4','chr12','chr13','chr19'),
                ID %in% vars) %>%
  dplyr::mutate(`#Chr` = c('chr1_variant1','chr1_variant2','chr4','chr12','chr13','chr19'))

refs = c('C','G','G','G','C','A') ; alts = c('T','A','A','A','T','G')
neas = c('T','A','A','G','C','A') ; eas = c('C','G','G','A','T','G')
gts = getAllGenotypes(refs,alts) #GRCh38 
color.set.1 = c('#fee8c8','#fdbb84','#e34a33')
color.set.2 = c('#deebf7','#9ecae1','#3182bd')
color.set.by.loc = list(color.set.2,color.set.2,color.set.1,color.set.1,color.set.1,
                        color.set.1)

curr.index = 6 #pick index between 1-6 for each sQTL
# set.seed(1)
nea = refs[[curr.index]] ; ea = alts[[curr.index]]
curr.excised.all = transpose((quant.excised %>% dplyr::filter(`#Chr`==names(gts)[[curr.index]]))[,3:ncol(quant.excised)],keep.names='Sample') %>%
  dplyr::rename(QE=V1)

#now make common dataset. Doing manually as names a bit different
#going by curr.excised.all as simpler sample names, so less code needed
plot.df = data.frame(Sample=character(),GT=character(),QuantExcised=numeric())
for (row in 1:nrow(curr.excised.all)) {
  match = which(stringr::str_detect(gts[[curr.index]]$Sample,curr.excised.all$Sample[[row]]))
  if (length(match)>0) 
    plot.df %<>% dplyr::add_row(Sample=curr.excised.all$Sample[[row]],
                         GT=gts[[curr.index]]$GT[[match]],
                         QuantExcised=curr.excised.all$QE[[row]])
}

#plot results for quant excised figure
plot.df %<>% dplyr::mutate(GT = ifelse(GT == "0|0",glue('{nea}{nea}'),GT)) %>%
  dplyr::mutate(GT = ifelse(GT == "0|1" | GT == "1|0",glue('{nea}{ea}'),GT)) %>%
  dplyr::mutate(GT = ifelse(GT == "1|1",glue('{ea}{ea}'),GT))

# manage case where the EA is not the alternate allele
# Want to reverse letter order for heterozygote, and reorder the factor
if (ea != eas[[curr.index]]) {
  plot.df %<>% dplyr::mutate(GT=ifelse(str_detect(GT,glue('{nea}{ea}')),glue('{ea}{nea}'),GT))
  plot.df$GT = factor(plot.df$GT,levels=c(glue('{ea}{ea}'),glue('{ea}{nea}'),
                                          glue('{nea}{nea}')))
}else{
  plot.df$GT = factor(plot.df$GT,levels=c(glue('{nea}{nea}'),glue('{nea}{ea}'),
                                          glue('{ea}{ea}')))
}

gt.counts = as.numeric(table(plot.df$GT))
nea = neas[[curr.index]] ; ea = eas[[curr.index]]

# for testing nonrandom samples for sashimi plots - exploring DPP9
# tmp = plot.df %>% arrange(GT,QuantExcised)
# tmp = tmp[c(131:159,sample(160:261,29),c(262:290)),]
# tmp = tmp[c(1:29,sample(160:261,29),c(262:290)),]
  
plt = ggplot(plot.df,aes(x=GT,y=QuantExcised,fill=GT)) + geom_violin() +
  scale_fill_manual(values=color.set.by.loc[[curr.index]]) +
  geom_jitter() + theme_bw() + theme(text=element_text(size=12)) +
  ylab('Normalized Intron Excision Ratio') + geom_boxplot(width=0.2,fill='white') +
  stat_boxplot(geom='errorbar',width=0.3,color='black') +
  xlab('Genotype') + theme(text=element_text(size=16),legend.position='none') +
  scale_x_discrete(limit=c(glue('{nea}{nea}'),glue('{nea}{ea}'),glue('{ea}{ea}')),
    labels=c(glue('{nea}{nea} (n={gt.counts[[1]]})'),
             glue('{nea}{ea} (n={gt.counts[[2]]})'),
             glue('{ea}{ea} (n={gt.counts[[3]]})'))) +
  ylim(c(-3,3)) + ggtitle(plot.titles[[curr.index]])
print(plt)
print(table(plot.df$GT))

# save plot to make uniform size
ggsave(glue('violin_{curr.index}.png'),plt,width=6,height=4,dpi=400)

print(c(nea,ea))
p1 = t.test((plot.df %>% dplyr::filter(GT==glue('{nea}{nea}')))$QuantExcised,
       (plot.df %>% dplyr::filter(GT==glue('{nea}{ea}')))$QuantExcised)$p.value
p2 = t.test((plot.df %>% dplyr::filter(GT==glue('{nea}{ea}')))$QuantExcised,
       (plot.df %>% dplyr::filter(GT==glue('{ea}{ea}')))$QuantExcised)$p.value
p.adjust(c(p1,p2),'BH')

##############################################
#Next make Sashimi plots
#First, need to get 20 randomly-selected samples for each genotype
makeGGSashimi.BamList(plot.df,names(gts)[[curr.index]],F,all.samples=T) #Bam files GRCh37
# makeTrackplot.BamList(plot.df,names(gts)[[curr.index]])

# for working with ggsashimi:
# go to line ~923 of ggsashimi.py to select the appropriate intron to highlight
# go to line ~1027 to put in the right order

################################################
# Get CPM by genotype
samtools.reads = getAllSamtoolsReads()

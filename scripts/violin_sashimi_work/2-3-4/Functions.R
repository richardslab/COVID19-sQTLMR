getAllGenotypes = function(neas,eas) {
  #transpose function from data.table
  
  chr = c('chr1_variant1','chr1_variant2','chr4','chr12','chr13','chr19')
  nea = neas
  ea = eas
  
  gt = list()
  
  for (c in chr) {
    print(glue('{c}_var.vcf.gz'))
    tmp = read.vcfR(glue('{c}_var.vcf.gz'))
    gt[[length(gt)+1]] = transpose(as.data.frame(extract.gt(tmp)) %>% dplyr::select_if(~ !any(is.na(.))),keep.names='Sample') %>%
      dplyr::rename(GT=V1)
  }
  names(gt) = chr
  return(gt)
}
getRandomSamplesForSashimi = function(gts) { #input is list of genotypes for each chromosome
  set.seed(1)
  samples = character()
  for (g in gts) {
    num.cases = ifelse(min(table(g$GT)) < 20,min(table(g$GT),20),20)
    samples %<>% append(c(sample((g %>% dplyr::filter(GT=='0|0'))$Sample,num.cases),
                          sample((g %>% dplyr::filter(GT=='0|1' | GT=='1|0'))$Sample,num.cases),
                          sample((g %>% dplyr::filter(GT=='1|1'))$Sample,num.cases)))
  }
  out.df = data.frame(Samples=samples) %>% dplyr::mutate(Samples = stringr::str_replace(Samples, "^(([^-]+-){1}[^-]+)-.*", "\\1"))
  return(out.df %>% dplyr::distinct())
}
runSamtoolsFunctions = function() {
  #first get all bam files
  bam.files = list.files('/project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung',
                         recursive=T,pattern='.bam',full.names=T,include.dirs=T)
  
  #next cut out files already processed or temp
  already.indexed.ids = bam.files[which(str_detect(bam.files,'sorted'))] %>%
    gsub('.*/','',.) %>% gsub('.sorted.*','',.) %>% unique()
  for (id in already.indexed.ids) bam.files = bam.files[which(!str_detect(bam.files,id))]
  bam.files = bam.files[!str_detect(bam.files,'temp')]
 
  #make necessary commands
  new.names = str_replace(bam.files,'.bam','.sorted.bam')
  commands = glue('echo `date` ; samtools sort {bam.files} -o {new.names} -@ 5 ; samtools index {new.names} ; echo `date`')
  
  write.table(commands,'samtools_sort_index_new.sh',row.names=F,col.names=F,quote=F)
}
makeGGSashimi.BamList = function(df,name,pick.extremes) {
  #get dictionary
  dictionary = vroom('/project/richards/restricted/dbGap/prj_32756/Documents/linking_file_GTEx.txt') %>%
    dplyr::filter(submitted_subject_id %in% df$Sample,body_site == 'Lung')
  bam.files = list.files('/project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung',
                         recursive=T,pattern='.sorted.bam',full.names=T,include.dirs=T)
  bam.files = bam.files[which(!str_detect(bam.files,'bai') & !str_detect(bam.files,'tmp'))]
  out.df = data.frame(ID=NA,Path=bam.files,Label=NA,Excis=NA) %>% 
    dplyr::filter(!str_detect(Path,'temp')) %>%
    dplyr::mutate(ID=gsub('.*/','',Path),ID=gsub('.sorted.bam','',ID))
  for (row in 1:nrow(out.df)) {
    dict.match = dictionary$submitted_subject_id[which(dictionary$Run == out.df$ID[[row]])]
    if (length(dict.match)>0) out.df$Label[[row]] = dict.match
    gt.match = as.character(df$GT[which(df$Sample == out.df$Label[[row]])])
    curr.excis = as.numeric(df$QuantExcised[which(df$Sample == out.df$Label[[row]])])
    if (length(gt.match)>0) {
      out.df$Label[[row]] = gt.match
      out.df$Excis[[row]] = curr.excis
    }else out.df$Label[[row]] = NA
  }
  out.df %<>% tidyr::drop_na(Label)
  
  #now to get same number of samples per genotype
  min.size = min(table(out.df$Label))
  indices.to.keep = numeric()
  for (g in unique(out.df$Label)) { #go by gt
    matching.indices = which(out.df$Label == g)
    indices.to.keep %<>% append(sample(matching.indices,min.size))
  }
  
  write.df = out.df[indices.to.keep,]
  
  vroom_write(write.df,glue('ggsashimi/{name}_input_bams_gtex.tsv'),col_names=F)
}

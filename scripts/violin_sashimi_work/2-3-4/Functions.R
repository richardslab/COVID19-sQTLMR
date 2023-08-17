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
makeGGSashimi.BamList = function(df,name,pick.extremes,all.samples) {
  #get dictionary
  dictionary = vroom('/project/richards/restricted/dbGap/prj_32756/Documents/linking_file_GTEx.txt') %>%
    dplyr::filter(submitted_subject_id %in% df$Sample,body_site == 'Lung')
  bam.files = list.files('Lung_GTEx_v8_Bam',
                         recursive=T,pattern='.bam',full.names=T,include.dirs=T)
  bam.files = bam.files[which(!str_detect(bam.files,'bai') & !str_detect(bam.files,'tmp'))]
  out.df = data.frame(ID=NA,Path=bam.files,Label=NA,Excis=NA) %>% 
    dplyr::filter(!str_detect(Path,'temp'),!str_detect(Path,'GTEX-111VG-0726-SM-5GIDC')) %>% #truncated file causing problems
    dplyr::mutate(ID=gsub('.*/','',Path),ID=gsub('.Aligned.sortedByCoord.out.patched.md.bam','',ID),
                  Path=glue('/scratch/richards/julian.willett/14.GTEx_sQTL_Tomoko_Work/{Path}'))
  for (row in 1:nrow(out.df)) {
    # dict.match = dictionary$submitted_subject_id[which(dictionary$biospecimen_repository_sample_id == out.df$ID[[row]])]
    # if (length(dict.match)>0) out.df$Label[[row]] = dict.match
    second.hyphen = str_locate_all(out.df$ID[[row]], "-")[[1]][2,'start']
    row.id = str_sub(out.df$ID[[row]], 1, second.hyphen-1)
    gt.match = as.character(df$GT[which(df$Sample == row.id)])
    curr.excis = as.numeric(df$QuantExcised[which(df$Sample == row.id)])
    if (length(gt.match)>0) {
      out.df$Label[[row]] = gt.match
      out.df$Excis[[row]] = curr.excis
    }else out.df$Label[[row]] = NA
  }
  print(names(out.df))
  out.df = out.df[which(!is.na(out.df$Label)),]

  #now to get same number of samples per genotype
  # min.size = 20
  # indices.to.keep = numeric()
  # for (g in unique(out.df$Label)) { #go by gt
  #   matching.indices = which(out.df$Label == g)
  #     
  #   if (!all.samples)
  #     indices.to.keep %<>% append(sample(matching.indices,min.size))
  #   else
  #     indices.to.keep %<>% append(matching.indices)
  # }
  # 
  # write.df = out.df[indices.to.keep,]
  write.df = out.df
  vroom_write(write.df,glue('ggsashimi/{name}_input_bams_gtex.tsv'),col_names=F)
  
  # next make command for getting local coverage at site:
  unique.gt = unique(write.df$Label)
  if (name == 'chr1_variant1') region = "chr1:155192209-155192809"
  else if (name == 'chr1_variant2') region = "chr1:156233609-156236509"
  else if (name == 'chr4') region = "chr4:105896843-105928843"
  else if (name == 'chr12') region = "chr12:112917600-112919595"
  else if (name == 'chr13') region = "chr13:112875786-112881876"
  else if (name == 'chr19') region = "chr19:4714037-4719988"
  
  for (label in unique.gt) { 
    gt.paths = (write.df %>% dplyr::filter(Label == label))$Path
    cmd = glue("samtools coverage -r {region} {paste(gt.paths,collapse = ' ')} > {name}_{label}_coverage.txt")
    print(cmd)
    writeLines(cmd,glue('coverage_command_{name}_{label}.sh'))
  }
}
makeTrackplot.BamList = function(df,name) {
  bam.files = list.files('Lung_GTEx_v8_Bam',
                         recursive=T,pattern='.sorted.bam',full.names=T,include.dirs=T)
  bam.files = bam.files[which(!str_detect(bam.files,'bai') & !str_detect(bam.files,'tmp'))]
  out.df = data.frame(SampleName=NA,Path=bam.files,FileType='bam',Label=NA,Color=NA,
                      Agg='fru') %>%
    mutate(SampleName=str_remove(Path, "\\..*"),SampleName=str_remove(SampleName,".*/")) %>%
    filter(!str_detect(Path,'temp'))
  dictionary = vroom('/project/richards/restricted/dbGap/prj_32756/Documents/linking_file_GTEx.txt') %>%
    dplyr::filter(submitted_subject_id %in% df$Sample,body_site == 'Lung') %>%
    mutate(SampleName=biospecimen_repository_sample_id)
  out.df = merge(out.df,dictionary)

  # add num reads
  num_reads = vroom('NumReads/bam_counts_region.txt',col_names=F) %>%
    mutate(X1 = str_remove(X1, ".*\\/"), X1 = str_remove(X1,'\\..*'))
  names(num_reads)[1:3] = c('SampleName','Region','NumReads')
  out.df.with.reads = merge(out.df,num_reads,by='SampleName') %>%
    select(SampleName,Region,Path,FileType,Label,Color,Agg,NumReads,SampleName)
  
  if (name == 'chr1_variant1') out.df.with.reads %<>% filter(Region=='muc1')
  else if (name == 'chr1_variant2') out.df.with.reads %<>% filter(Region=='pmf1')
  else if (name == 'chr4') out.df.with.reads %<>% filter(Region=='npnt')
  else if (name == 'chr12') out.df.with.reads %<>% filter(Region=='oas1')
  else if (name == 'chr13') out.df.with.reads %<>% filter(Region=='atp11a')
  else if (name == 'chr19') out.df.with.reads %<>% filter(Region=='dpp9')
  out.df.with.reads %<>% distinct()

  for (row in 1:nrow(out.df.with.reads)) {
    if (name == 'chr1_variant1') colors = c('#3182bd','#9ecae1','#deebf7')
    else if (name == 'chr1_variant2') colors = c('#3182bd','#deebf7','#9ecae1')
    else if (name == 'chr4') colors = c('#e34a33','#fdbb84','#fee8c8')
    else if (name == 'chr12') colors = c('#e34a33','#fdbb84','#fee8c8')
    else if (name == 'chr13') colors = c('#fdbb84','#fee8c8','#e34a33')
    else if (name == 'chr19') colors = c('#fdbb84','#fee8c8','#e34a33')
    genotype.vector = c(paste0(nea,nea),paste0(nea,ea),paste0(ea,ea))

    cut.id = unlist(strsplit(out.df.with.reads$SampleName[[row]], "-"))
    just.id = paste0(cut.id[1], "-", cut.id[2])
    gt.match = as.character(df$GT[which(df$Sample == just.id)])
    if (length(gt.match)>0) {
      out.df.with.reads$Label[[row]] = gt.match
      out.df.with.reads$Color[[row]] = colors[which(gt.match == genotype.vector)]
    }else out.df$Label[[row]] = NA
  }
  write.df = out.df.with.reads %>% filter(!is.na(Label),!is.na(NumReads)) %>%
    select(-SampleName,-Region)

  vroom_write(write.df,glue('trackplot/{name}_input_bams_gtex.tsv'),col_names=F)
  write.df
}
getAllSamtoolsReads = function() {
  files = list.files('deepTools/reads/',full.names=T,include.dirs=T)
  df = data.frame(Filepath=character(),ID=character(),Chr1_Var1_Reads=numeric(),
                  Chr1_Var2_Reads=numeric(),Chr4_Reads=numeric(),Chr12_Reads=numeric(),
                  Chr13_Reads=numeric(),Chr19_Reads=numeric(),TotalReads=numeric())
  
  # go through files and get data
  for (f in files) {
    tmp = vroom(f,skip=1,col_names=F,show_col_types=F)
    names(tmp)= c('Filepath','ID','Chr1_Var1_Reads','Chr1_Var2_Reads','Chr4_Reads',
                  'Chr12_Reads','Chr13_Reads','Chr19_Reads','TotalReads')
    df %<>% add_row(tmp)
  }
  
  # calculate CPM
  df %<>% mutate(Chr1_Var1_CPM=Chr1_Var1_Reads/TotalReads*1e6,
                 Chr1_Var2_CPM=Chr1_Var2_Reads/TotalReads*1e6,
                 Chr4_CPM=Chr4_Reads/TotalReads*1e6,
                 Chr12_CPM=Chr12_Reads/TotalReads*1e6,
                 Chr13_CPM=Chr13_Reads/TotalReads*1e6,
                 Chr19_CPM=Chr19_Reads/TotalReads*1e6,
                 ID=str_replace(ID,'\\.',''))
  
  # add genotype information
  df %<>% mutate(Chr1_Var1_GT=NA,Chr1_Var2_GT=NA,Chr4_GT=NA,Chr12_GT=NA,Chr13_GT=NA,
                 Chr19_GT=NA) %>% filter(!is.na(TotalReads))
  chr1.var1.file = vroom('ggsashimi/chr1_variant1_input_bams_gtex.tsv',col_names=F)
  chr1.var2.file = vroom('ggsashimi/chr1_variant2_input_bams_gtex.tsv',col_names=F)
  chr4.file = vroom('ggsashimi/chr4_input_bams_gtex.tsv',col_names=F)
  chr12.file = vroom('ggsashimi/chr12_input_bams_gtex.tsv',col_names=F)
  chr13.file = vroom('ggsashimi/chr13_input_bams_gtex.tsv',col_names=F)
  chr19.file = vroom('ggsashimi/chr19_input_bams_gtex.tsv',col_names=F)
  
  for (r in 1:nrow(df)) {
    print(r)
    df$Chr1_Var1_GT[[r]] = try(chr1.var1.file[[which(chr1.var1.file$X1 == df$ID[[r]]),'X3']])
    df$Chr1_Var2_GT[[r]] = try(chr1.var2.file[[which(chr1.var2.file$X1 == df$ID[[r]]),'X3']])
    df$Chr4_GT[[r]] = try(chr4.file[[which(chr4.file$X1 == df$ID[[r]]),'X3']])
    df$Chr12_GT[[r]] = try(chr12.file[[which(chr12.file$X1 == df$ID[[r]]),'X3']])
    df$Chr13_GT[[r]] = try(chr13.file[[which(chr13.file$X1 == df$ID[[r]]),'X3']])
    df$Chr19_GT[[r]] = try(chr19.file[[which(chr19.file$X1 == df$ID[[r]]),'X3']])
  }
  df %<>% filter(!str_detect(Chr1_Var1_GT,'Error'))
  
  # print out the average CPM by locus and genotype
  for (gt in unique(df$Chr1_Var1_GT)) {
    tmp.df = df %>% filter(Chr1_Var1_GT == gt)
    print(glue('Chr1_Var1 GT: {gt}. N: {nrow(tmp.df)}. Mean CPM: {mean(tmp.df$Chr1_Var1_CPM)}'))
  }
  for (gt in unique(df$Chr1_Var2_GT)) {
    tmp.df = df %>% filter(Chr1_Var2_GT == gt)
    print(glue('Chr1_Var2 GT: {gt}. N: {nrow(tmp.df)}. Mean CPM: {mean(tmp.df$Chr1_Var2_CPM)}'))
  }
  for (gt in unique(df$Chr4_GT)) {
    tmp.df = df %>% filter(Chr4_GT == gt)
    print(glue('Chr4 GT: {gt}. N: {nrow(tmp.df)}. Mean CPM: {mean(tmp.df$Chr4_CPM)}'))
  }
  for (gt in unique(df$Chr12_GT)) {
    tmp.df = df %>% filter(Chr12_GT == gt)
    print(glue('Chr12 GT: {gt}. N: {nrow(tmp.df)}. Mean CPM: {mean(tmp.df$Chr12_CPM)}'))
  }
  for (gt in unique(df$Chr13_GT)) {
    tmp.df = df %>% filter(Chr13_GT == gt)
    print(glue('Chr13 GT: {gt}. N: {nrow(tmp.df)}. Mean CPM: {mean(tmp.df$Chr13_CPM)}'))
  }
  for (gt in unique(df$Chr19_GT)) {
    tmp.df = df %>% filter(Chr19_GT == gt)
    print(glue('Chr19 GT: {gt}. N: {nrow(tmp.df)}. Mean CPM: {mean(tmp.df$Chr19_CPM)}'))
  }
  
  return(df)
}

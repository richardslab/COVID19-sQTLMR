import pandas as pd

exppath = "/project/richards/public/gtex/scratch/v8/gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_sQTL_all_associations"
outpath = "/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR/"

map = pd.read_csv(outpath+"sQTL_WBC.map", sep="\t")
for index, row in map.iterrows():
  dat = pd.read_parquet(exppath+"/Whole_Blood.v8.EUR.sqtl_allpairs."+row['CHR']+".parquet")
  select_dat = dat[dat['phenotype_id'].str.contains(row['group_id'])]
  select_dat.to_csv(outpath+"sQTL_WBC/"+row['group_id']+".tsv",sep="\t",index=False)

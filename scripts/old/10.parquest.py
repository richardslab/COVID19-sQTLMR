import pandas as pd

exppath = "/project/richards/public/gtex/scratch/v8/gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations"
outpath = "/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR/"

map = pd.read_csv(outpath+"eQTL.map", sep="\t")
for index, row in map.iterrows():
  dat = pd.read_parquet(exppath+"/Lung.v8.EUR.allpairs."+row['CHR']+".parquet")
  select_dat = dat[dat['phenotype_id'].str.contains(row['exposure'])]
  select_dat.to_csv(outpath+"eQTL_Lung/"+row['exposure']+".tsv",sep="\t",index=False)

import pandas as pd

exppath = "/project/richards/public/gtex/scratch/v8/gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_sQTL_all_associations"
outpath = "/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR/"

map = pd.read_csv(outpath+"HGI_r7_loci", sep="\t", header=None)
for index, row in map.iterrows():
  dat = pd.read_parquet(exppath+"/Whole_Blood.v8.EUR.sqtl_allpairs.chr"+str(row[0])+".parquet")
  dat1 = dat['variant_id'].str.split("_", expand=True)
  dat2 = pd.concat([dat, dat1], axis=1)
  dat2 = dat2.rename(columns={1: "POS"})
  dat2['POS'] = pd.to_numeric(dat2["POS"])
  select_dat = dat2[dat2['POS'] > int(row[1]) - 500000]
  select_dat = select_dat[select_dat['POS'] < int(row[1]) + 500000]
  select_dat.to_csv(outpath+"sQTL_WBC/"+str(row[2])+".tsv",sep="\t",index=False)

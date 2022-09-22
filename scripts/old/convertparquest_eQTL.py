import pandas as pd

exppath = "/project/richards/public/gtex/scratch/v8/gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/"
outpath = "/home/richards/tomoko.nakanishi/scratch/DATA/eQTL/GTEx"

for index in range(1, 23):
  dat = pd.read_parquet(exppath+"/Whole_Blood.v8.EUR.allpairs.chr"+str(index)+".parquet")
  dat.to_csv(outpath+"Whole_Blood.v8.EUR.allpairs.chr"+str(index)+".tsv.gz",sep="\t",index=False,compression='gzip')

for index in range(1, 23):
  dat = pd.read_parquet(exppath+"/Lung.v8.EUR.allpairs.chr"+str(index)+".parquet")
  dat.to_csv(outpath+"Lung.v8.EUR.allpairs.chr"+str(index)+".tsv.gz",sep="\t",index=False,compression='gzip')


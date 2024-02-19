# COVID19-sQTLMR

The codes and results repository for the MR analyses to understand the effect of alternative splicing on COVID-19.

**citation:** [Nakanishi T*†, Willet J*, _et al_. Alternative splicing in lung influences COVID-19 severity and respiratory diseases. _Nat Commun_ 2023](https://www.nature.com/articles/s41467-023-41912-4)

under `scripts`

* `01.makeoutcome.rds.R`: scripts to format outcome GWAS dataset to run MR.
* `02.makeexposure_sQTL_EUR.R`: scripts to format exposure GWAS dataset to run MR.
* `03.MRsQTL_EUR.R`: scripts to run MR using sQTL dataset as exposure. 
* `04.MRsQTL_summarise.R`: scripts to summarise MR results using sQTL dataset.
* `05.Fig3.MR.R`: scripts to make Fig 3 (forrest plot).
* `06.B1_MR.R`: scripts to run MR for B1 phenotype.
* `07.coloc.lung.R` & `08.coloc.wbc.R` : scripts to run colocalization analysis using sQTL dataset.
* `09.MRsQTL_coloc_sensitivity.R`: scripts to format colocalization analysis results.
* `10.MReQTL_EUR.R`: scripts to run MR using eQTL dataset as exposure.
* `11.coloc.wbc.eQTL.R` & `12.coloc.lung.eQTL.R`: scripts to run colocalization analysis using eQTL dataset.
* `13.MReQTL_sensitivity.R`: scripts to summarise MR results using eQTL dataset.
* `14.HPA.scRNAseq.heatmap.R`: scripts to make Fig. 4A and 4B.
* `15.scRNAseq_covid.R`: scripts to make Fig. 4C.
* `16.coloc_otherdiseases.R`: scripts to run colocalization for other diseases than COVID-19.

under `scripts/violin_sashimi_work`

Codes to generate Fig 2 (violin plot) and Supplementary Fig 1 (sashimi plot) were stored. For details, please refer to the Readme file within the folder.

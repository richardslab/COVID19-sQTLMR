#!/bin/bash

SCRATCH=/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR/genotype
IMPUTE_DIR=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/07.TopMed

for CHR in {1..22}
do
plink2 --bgen ${IMPUTE_DIR}/${CHR}.bgen 'ref-first' \
--sample ${IMPUTE_DIR}/1-23.sample \
--export A \
--extract ${SCRATCH}/leadsnp --out ${SCRATCH}/${CHR}
done

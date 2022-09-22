#!/bin/bash

WIRKDIR=/home/richards/tomoko.nakanishi/scratch/09.COVID19/12.COVID19-pQTLMR

mkdir -p ${WIRKDIR}/LD/

cat ${WIRKDIR}/sig_rsid | while read line
do
CHR=$(echo $line | cut -f2 -d " " | sed -e "s/chr//g")
rsid=$(echo $line | cut -f1 -d " " )
plink --bfile /project/richards/public/1000Genomes/scratch/wgs/phase3/plink/$CHR \
--keep /project/richards/public/1000Genomes/scratch/wgs/phase3/integrated_call_samples_v3.20130502.ALL.panel.EURonly \
--r2 \
--ld-snp $rsid \
--ld-window-kb 1000  \
--ld-window 99999 \
--ld-window-r2 0 \
--out $WIRKDIR/LD/$rsid
done

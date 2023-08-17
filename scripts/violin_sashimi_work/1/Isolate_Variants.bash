#!/bin/bash

#GRCh38

bcftools view \
GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
--regions chr1:155192276 --output-type z --output Output/chr1_variant1_var.vcf.gz

bcftools view \
GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
--regions chr1:156236330 --output-type z --output Output/chr1_variant2_var.vcf.gz

bcftools view \
GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
--regions chr4:105897896 --output-type z --output Output/chr4_var.vcf.gz

bcftools view \
GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
--regions chr12:112919388 --output-type z --output Output/chr12_var.vcf.gz

bcftools view \
GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
--regions chr13:112881427 --output-type z --output Output/chr13_var.vcf.gz

bcftools view \
GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
--regions chr19:4717660 --output-type z --output Output/chr19_var.vcf.gz

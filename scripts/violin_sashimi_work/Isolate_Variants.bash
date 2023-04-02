export PATH=~/projects/richards/share/bin:$PATH

bcftools view \
GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz \
--regions 1:155162067 --output-type z --output Output/chr1_variant1_var.vcf.gz

bcftools view \
GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz \
--regions 1:156206121 --output-type z --output Output/chr1_variant2_var.vcf.gz

bcftools view \
GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz \
--regions 4:106819053 --output-type z --output Output/chr4_var.vcf.gz

bcftools view \
GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz \
--regions 12:113357193 --output-type z --output Output/chr12_var.vcf.gz

bcftools view \
GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz \
--regions 13:113535741 --output-type z --output Output/chr13_var.vcf.gz

bcftools view \
GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz \
--regions 19:4717672 --output-type z --output Output/chr19_var.vcf.gz

#exome gatk output is grch37, not grch38
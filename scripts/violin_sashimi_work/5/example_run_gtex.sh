## ====================
##  ggsashimi examples
## ====================

#interpreting sashimi plots: https://miso.readthedocs.io/en/fastmiso/sashimi.html

module load StdEnv/2020 gcc/9.3.0 r/4.1.2 r-bundle-bioconductor/3.14 python/3.9.6

# GRCh38 coordinates, step up from grch37
bam_input=chr1_variant1_input_bams_gtex.tsv ; curr_coord="chr1:155192209-155192809" ; MC=1 ; palette=palette2.txt ; prfx="muc1" #MUC1
bam_input=chr1_variant2_input_bams_gtex.tsv ; curr_coord="chr1:156233609-156236509" ; MC=1 ; palette=palette1_2.txt ; prfx="pmf1" #PMF1
bam_input=chr4_input_bams_gtex.tsv ; curr_coord="chr4:105896843-105928843" ; MC=1 ; palette=palette4.txt ; prfx="npnt"
bam_input=chr12_input_bams_gtex.tsv ; curr_coord="chr12:112917600-112919595" ; MC=1 ; palette=palette12.txt ; prfx="oas1"
###bam_input=chr13_input_bams_gtex.tsv ; curr_coord="13:113530100-113535060" ; MC=10 ; palette=palette19.txt 
bam_input=chr13_input_bams_gtex.tsv ; curr_coord="chr13:112875786-112881876" ; MC=1 ; palette=palette13.txt ; prfx="atp11a"
###bam_input=chr19_input_bams_gtex.tsv ; curr_coord="chr19:4714049-4717700" ; MC=1 ; palette=palette19.txt #DPP9
bam_input=chr19_input_bams_gtex.tsv ; curr_coord="chr19:4714037-4719988" ; MC=1 ; palette=palette19.txt ; prfx="dpp9" #DPP9 extended


#go to line 1026 to adjust order

## Example #1. Overlay, intron shrinkage, gene annotation, PDF output, custom size and colors
echo `date`
../ggsashimi_${prfx}.py -o $prfx -b ../${bam_input} -c ${curr_coord} -A mean_j -g gencode.v43.chr_patch_hapl_scaff.annotation.gtf -M $MC -C 3 -O 3 --shrink --alpha 0.25 --base-size=32 --ann-height=6 --height=3 --width=18 -P $palette
echo $prfx `date`

# for simple annotation (for oas1)
../ggsashimi_${prfx}.py -o $prfx -b ../${bam_input} -c ${curr_coord} -A mean_j -g gencode.v43.basic.annotation.gtf -M $MC -C 3 -O 3 --shrink --alpha 0.25 --base-size=32 --ann-height=6 --height=3 --width=18 -P $palette


## Example #2. Mean coverage and number of reads supporting inclusion and exclusion, no gene annotation, TIFF output (350 PPI), custom size, default colors, fixed y-scale
#../ggsashimi.py -b ../${bam_input} -c ${curr_coord} -M $MC -C 3 -O 3 -A mean --alpha 1 -F tiff -R 700 --base-size=28 --height=3 --width=18 --fix-y-scale


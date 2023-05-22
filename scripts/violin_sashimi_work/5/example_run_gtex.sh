## ====================
##  ggsashimi examples
## ====================

#interpreting sashimi plots: https://miso.readthedocs.io/en/fastmiso/sashimi.html

export PATH=~/miniconda3/bin:$PATH
source activate ggsashimi
module load r/4.1.2 r-bundle-bioconductor/3.14

bam_input=chr1_variant1_input_bams_gtex.tsv ; curr_coord="1:155162000-155162600" ; MC=30 ; palette=palette2.txt #MUC1
bam_input=chr1_variant2_input_bams_gtex.tsv ; curr_coord="1:156203400-156206300" ; MC=10 ; palette=palette1_2.txt #PMF1
bam_input=chr4_input_bams_gtex.tsv ; curr_coord="4:106818000-106850000" ; MC=10 ; palette=palette12.txt
bam_input=chr12_input_bams_gtex.tsv ; curr_coord="12:113355405-113357400" ; MC=1 ; palette=palette12.txt
###bam_input=chr13_input_bams_gtex.tsv ; curr_coord="13:113530100-113535060" ; MC=10 ; palette=palette19.txt
bam_input=chr13_input_bams_gtex.tsv ; curr_coord="13:113530100-113536190" ; MC=10 ; palette=palette19.txt
###bam_input=chr19_input_bams_gtex.tsv ; curr_coord="19:4714049-4718227" ; MC=3 ; palette=palette19.txt #DPP9
bam_input=chr19_input_bams_gtex.tsv ; curr_coord="19:4714049-4720000" ; MC=3 ; palette=palette19.txt #DPP9 extended

#go to line 1026 to adjust order

## Example #1. Overlay, intron shrinkage, gene annotation, PDF output, custom size and colors
../ggsashimi.py -b ../${bam_input} -c ${curr_coord} -A mean_j -g Homo_sapiens.GRCh37.87.gtf -M $MC -C 3 -O 3 --shrink --alpha 0.25 --base-size=32 --ann-height=6 --height=3 --width=18 -P $palette
## Example #2. Mean coverage and number of reads supporting inclusion and exclusion, no gene annotation, TIFF output (350 PPI), custom size, default colors, fixed y-scale
#../ggsashimi.py -b ../${bam_input} -c ${curr_coord} -M $MC -C 3 -O 3 -A mean --alpha 1 -F tiff -R 700 --base-size=28 --height=3 --width=18 --fix-y-scale



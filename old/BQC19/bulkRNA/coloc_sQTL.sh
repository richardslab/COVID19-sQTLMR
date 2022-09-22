#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=rrg-brichard
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=tomoko.nakanishi@mail.mcgill.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL
#SBATCH --array=1-49


module load r/4.1.2 

for PHENO in A2 B2 C2
do
CHR=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/tomoco/scratch/coloc/HGI_r7_loci | awk '{print $1}')
POS=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/tomoco/scratch/coloc/HGI_r7_loci | awk '{print $2}')
LOC=$(sed -n ${SLURM_ARRAY_TASK_ID}p /home/tomoco/scratch/coloc/HGI_r7_loci | awk '{print $3}')
Rscript /home/tomoco/projects/def-brichard/tomoco/COVID19/coloc/01.coloc_sQTL.R ${PHENO} ${CHR} ${POS} ${LOC}
done

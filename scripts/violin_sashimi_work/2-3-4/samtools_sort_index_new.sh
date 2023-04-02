#!/bin/bash
#PBS -N samtools_sqtl
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=16G
#PBS -l vmem=16G

export PATH=~/miniconda3/bin:$PATH
source activate ggsashimi

echo `date` ; samtools sort /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1446546/SRR1446546.bam -o /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1446546/SRR1446546.sorted.bam -@ 5 ; samtools index /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1446546/SRR1446546.sorted.bam ; echo `date`
echo `date` ; samtools sort /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1454415/SRR1454415.bam -o /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1454415/SRR1454415.sorted.bam -@ 5 ; samtools index /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1454415/SRR1454415.sorted.bam ; echo `date`
echo `date` ; samtools sort /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1471050/SRR1471050.bam -o /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1471050/SRR1471050.sorted.bam -@ 5 ; samtools index /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1471050/SRR1471050.sorted.bam ; echo `date`
echo `date` ; samtools sort /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1472462/SRR1472462.bam -o /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1472462/SRR1472462.sorted.bam -@ 5 ; samtools index /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1472462/SRR1472462.sorted.bam ; echo `date`
echo `date` ; samtools sort /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1478015/SRR1478015.bam -o /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1478015/SRR1478015.sorted.bam -@ 5 ; samtools index /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR1478015/SRR1478015.sorted.bam ; echo `date`
echo `date` ; samtools sort /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR2184865/SRR2184865.bam -o /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR2184865/SRR2184865.sorted.bam -@ 5 ; samtools index /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR2184865/SRR2184865.sorted.bam ; echo `date`
echo `date` ; samtools sort /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR2185097/SRR2185097.bam -o /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR2185097/SRR2185097.sorted.bam -@ 5 ; samtools index /project/richards/restricted/dbGap/prj_32756/RNA_Seq/Lung/SRR2185097/SRR2185097.sorted.bam ; echo `date`

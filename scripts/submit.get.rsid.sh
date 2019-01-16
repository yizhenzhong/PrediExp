#!/bin/bash
#MSUB -A b1042
#MSUB -l walltime=4:00:00
#MSUB -l nodes=1:ppn=20
#MSUB -j oe
#MSUB -q genomics
#MSUB -N dexseq_peer
cd /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/scripts
module load shapeit
module load plink
module load vcftools
module load python
#module load R
module load  R/3.4.4
module load java
python get_rsid.py \
        /projects/b1047/zhong/software/utility/rsid_SNP_b150.txt \
        /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/annotation3.txt \
        /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/ annotation_snp

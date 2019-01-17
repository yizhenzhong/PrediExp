#!/bin/bash
#MSUB -A b1042
#MSUB -l walltime=48:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -j oe
#MSUB -q genomics
#MSUB -N training

#cd /projects/b1047/zhong/v7_gtex/genotype/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1
#cd /projects/b1047/zhong/Hepatocyte_project/script
#cd /projects/b1047/zhong/Hepatocyte_project/results/multiQTL
cd /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediExp/PredictDB_Pipeline_GTEx_v7/model_training/scripts

module load shapeit
module load plink
module load vcftools
module load python
module load R

Rscript AA_tiss_chrom_training_new_genotype.R $MOAB_JOBARRAYINDEX


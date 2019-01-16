#!/bin/bash
#MSUB -A b1042
#MSUB -l walltime=48:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -j oe
#MSUB -q genomics
#MSUB -N peerexon

#cd /projects/b1047/zhong/v7_gtex/genotype/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1
#cd /projects/b1047/zhong/Hepatocyte_project/script
#cd /projects/b1047/zhong/Hepatocyte_project/results/multiQTL
cd /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/PredictDB_Pipeline_GTEx_v7/model_training/scripts

module load shapeit
module load plink
module load vcftools
module load python
module load R

#Rscript gtex_tiss_chrom_training.R $MOAB_JOBARRAYINDEX

#Rscript run_PEER.r /projects/b1047/zhong/Hepatocyte_project/Yilin/baseline_hg19_exon.bed.gz baseline_hg19_exon 15 --output_dir /projects/b1047/zhong/Hepatocyte_project/Yilin/

#Rscript /projects/b1047/zhong/Hepatocyte_project/script/Step4CalleQTL.R
##Rscript Step3ModelParameter.R
##python subset_rows.py 
#Rscript run_PEER_withCov.r ../data/expression/condition1_expression_tmm_normalization.csv ../data/expression/condition1_expression_tmm_norm_cov.csv condition1_expression_tmm_normalization_withcov 15 -o ../data/expression/PEER/ -n $MOAB_JOBARRAYINDEX
#Rscript countclust.R --k $MOAB_JOBARRAYINDEX
#/projects/b1047/zhong/Hepatocyte_project/script/pipeline.sh 

#python ../../../software/ancestry_pipeline-master/get_BY_BH_correction.py >linear.log
#Rscript run_PEER.r ../data/expression/condition1_expression_tmm_normalization.csv  condition1_expression_tmm_normalization 15 -o ../data/expression/PEER/ -n $MOAB_JOBARRAYINDEX

#Rscript ../../software/ancestry_pipeline-master/eqtl.matrixEL.r
#Rscript ../../software/ancestry_pipeline-master/eqtl.matrixEL.r ../data/genotype/eQTL/condition1_expression_tmm_normalization_dosage_v7.txt ../data/genotype/eQTL/hepatocyte_snp.position ../data/genotype/eQTL/condition1_expression_tmm_normalization_exp_v7.txt ../data/genotype/eQTL/Hepatocyte_gene_pos.txt  ../data/genotype/eQTL/condition1_10peer_covariate.txt ../data/genotype/eQTL/output/condition1_10peer_covariate.txt
Rscript gtex_tiss_chrom_training.R 4

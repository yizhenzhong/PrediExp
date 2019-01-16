#!/bin/bash
#MSUB -A b1042
#MSUB -l walltime=4:00:00
#MSUB -l nodes=1:ppn=20
#MSUB -j oe
#MSUB -q genomics
#MSUB -N dexseq_peer
cd /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/scripts
#cd /projects/b1047/zhong/Hepatocyte_project/data/genotype
#cd /projects/b1047/zhong/v7_gtex/genotype/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1
#cd /projects/b1047/zhong/Hepatocyte_project/script
#cd /projects/b1047/zhong/Hepatocyte_project/results/multiQTL
#cd /projects/b1047/zhong/Hepatocyte_project/meta-tissue/data
module load shapeit
module load plink
module load vcftools
module load python
#module load R
module load  R/3.4.4
module load java
#module load blas-lapack/3.6.0_gcc
#/projects/b1047/zhong/software/EIG-7.2.1/bin/convertf -p par.PED.EIGENSTRAT
python get_rsid.py /projects/b1047/zhong/software/utility/rsid_SNP_b150.txt /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/annotation3.txt /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/ annotation4
#Rscript run_PEER.r /projects/b1047/zhong/Hepatocyte_project/Yilin/Exon/dexseq_baseline_hg19_exon.bed.gz dexseq_baseline_hg19_exon 15 --output_dir /projects/b1047/zhong/Hepatocyte_project/Yilin/Exon/
#./run_step3.sh
#Rscript phi1_replication.R
#Rscript Step3ModelParameter.R##Rscript Step4CalleQTL.R
#script combineResults.R
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

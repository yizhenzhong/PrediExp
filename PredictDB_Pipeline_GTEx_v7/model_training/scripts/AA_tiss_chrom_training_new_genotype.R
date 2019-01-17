setwd("/projects/b1047/zhong/Hepatocyte_project/Yilin/PrediExp/PredictDB_Pipeline_GTEx_v7/model_training/scripts")
source("gtex_v7_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
snp_annot_file <- "../../../data/input/annotations/snp_annotation/annotation_snp.chr"%&% chrom %&% ".txt"
gene_annot_file <- "../../../data/intermediate/annotations/gene_annotation/gencode.v25lift37.annotation.updateid.txt"
genotype_file <- "../../../data/input/genotypes/genotype_filter.chr"%&% chrom %&% ".txt"
expression_file <- "../../../data/intermediate/expression_phenotypes/condition1_exp.txt"
covariates_file <- "../../../data/intermediate/expression_phenotypes/condition1_cov.txt"
prefix <- "filtered_nested_cv"
print(chrom)
print(snp_annot_file)
main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)



for i in {1..22}; do
        python split_genotype_by_chr.py /projects/b1047/zhong/Hepatocyte_project/data/genotype/eQTL/condition1_expression_tmm_normalization_dosage_v7_2.txt  \
                /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/genotypes/genotype  \
                /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/snp_annotation.chr$i.txt;
done        

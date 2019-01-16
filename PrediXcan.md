
## PrediXcan

### To get the SNP annotation file
To provide SNPID, pos, ref, and alt allele information

- add rsid from dbsnp150
Note: 
*37114456 SNPs in dbsnp150*
*425583 out of 7034069 SNPs with no rsIDs/bialleic SNPs, and these SNPs will be ignored during the training

https://s3.amazonaws.com/predictdb2/GTEx-V7_HapMap-2017-11-29_README.txt
*SNPs were then filtered to those
with minor allele frequency >= 1%, R2 >= 0.8, biallelic only, and those which
would be unambiguously stranded, i.e. the polymorphism is not to the base's
complement.*
```bash
python get_rsid.py \
/projects/b1047/zhong/software/utility/rsid_SNP_b150.txt \
/projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/annotation3.txt
/projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/ \
annotation_snp


[yzu280@quser11 snp_annotation]$ awk '$6=="NA" {print}' annotation_snp.txt | wc -l
425583

```
- make sure the reference and alternative allele matches with the dosage files.
[yzu280@quser13 scripts]$ python check_allele.py

### split the SNP files by chr
```python
#Skip non-single letter polymorphisms
 53             if len(refAllele) > 1 or len(effectAllele) > 1:
 54                 continue
 55             # Skip ambiguous strands
 56             if SNP_COMPLEMENT[refAllele] == effectAllele:
 57                 continue
 58             if rsid == 'NA':
 59                 continue
```
```bash
[yzu280@quser11 scripts]$ python split_snp_annot_by_chr.py  \
/projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/annotation_snp.txt  \
/projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/annotation_snp
```
Before filtering, the number of SNPs is 7034069
After filtering, the number of SNPs is 5189842


### split the genotype by chr
Note: need to use the annotation files after filtering
```bash
python split_genotype_by_chr.py /projects/b1047/zhong/Hepatocyte_project/data/genotype/eQTL/condition1_expression_tmm_normalization_dosage_v7_2.txt \ /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/genotypes/genotype_filter \
/projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/input/annotations/snp_annotation/annotation_snp.filterd.txt
```

 - [x] the genotype number is 5189842

###  get the gene expression annotations/covariate file
- Transpose the expression and covariates
- Covariates include: 10 PEER, sex, platform, batch, PC1-3


### Start the training





<!---


#### Test the PrediXcan model built with white people
```bash
python convert_plink_to_dosage.py \
-b /projects/b1047/zhong/Hepatocyte_project/data/genotype/baseline.genotype.eqtl.excludeoutlier.header \
-o /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/aa_input/aa_genotype -p plink


python PrediXcan.py --predict \
--dosages /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/aa_input/ \
--dosages_prefix aa_genotype \
--sample /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/aa_input/aa_genotype.nosex  \
--weights ../Liver/gtex_v7_Liver_imputed_europeans_tw_0.5_signif.db \
--output_prefix /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/aa_input/

```



```R
> quantile(cor_exp, c(.0,.25,.5,.75,1), na.rm = T)
        0%        25%        50%        75%       100% 
-0.3752232 -0.0263705  0.0792258  0.1993307  0.9068370 
# > length(exp_id)
# [1] 23042
# > length(predicted_id)
# [1] 3351
```
Note: The correlation is low.
possible explanations are:
1. less SNPs are used for GTEx model
2. Should compare the expression after regressing out the covariates.
3. 
after regressing out the covariates, the correlation is actually less
```r
> quantile(cor_exp, c(.0,.25,.5,.75,1), na.rm = T)
         0%         25%         50%         75%        100% 
-0.38839424 -0.02202818  0.08649508  0.22052888  0.76986455 
#
```
### Test the PrediXcan model built on white people on GTEx Liver AA people


```bash
python convert_plink_to_dosage.py \
-b /projects/b1047/zhong/v7_gtex/genotype/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_maf0.05_liver \
-o /projects/b1047/zhong/Hepatocyte_project/Yilin/PrediXcan/data/aa_input/gtex_liver \
-p plink
```


-->





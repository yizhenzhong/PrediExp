for chrom in {1..22}
do
        echo -e "    $chrom";
        Rscript ../../scripts/genotypes_to_transposed_RDS.R \
        ../../data/input/genotypes/condition1_expression.chr$chrom.txt \
        ../../data/intermediate/genotypes/condition1_expression.chr$chrom;
        sleep .5;
done

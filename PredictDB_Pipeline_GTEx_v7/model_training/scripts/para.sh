for j in {1..22}; do
        echo $j
        cp submit.sh submit_rpkm_$j.sh
        echo "Rscript gtex_tiss_chrom_training_rpkm.R $j" >> submit_rpkm_$j.sh        
        msub submit_rpkm_$j.sh

        #Rscript gtex_tiss_chrom_training_rpkm.R $j
done


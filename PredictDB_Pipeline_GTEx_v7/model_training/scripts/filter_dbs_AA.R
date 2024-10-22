library(dplyr)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b,sep='')
argv <- commandArgs(trailingOnly=T)
prefix <- argv[1]
driver <- dbDriver("SQLite")

unfiltered_db <-'../dbs/' %&% prefix %&% '.db' 
filtered_db <- '../dbs/' %&% prefix %&% 'filtered.db' 
in_conn <- dbConnect(driver, unfiltered_db)
out_conn <- dbConnect(driver, filtered_db)
model_summaries <- dbGetQuery(in_conn, 'select * from model_summaries where zscore_pval < 0.05 and rho_avg > 0.1')
model_summaries <- model_summaries %>% dplyr::rename(pred.perf.R2 = rho_avg_squared, genename = gene_name, pred.perf.pval = zscore_pval, n.snps.in.model = n_snps_in_model)
print(dim(model_summaries))
model_summaries$pred.perf.qval <- NA
dbWriteTable(out_conn, 'extra', model_summaries, overwrite=TRUE)
construction <- dbGetQuery(in_conn, 'select * from construction')
dbWriteTable(out_conn, 'construction', construction, overwrite=TRUE)
print(dbReadTable(out_conn, "construction"))
sample_info <- dbGetQuery(in_conn, 'select * from sample_info')
dbWriteTable(out_conn, 'sample_info', sample_info, overwrite=TRUE)
weights <- dbGetQuery(in_conn, 'select * from weights')
weights <- weights %>% dplyr::filter(gene %in% model_summaries$gene) %>% dplyr::rename(eff_allele = alt, ref_allele = ref, weight = beta)
dbWriteTable(out_conn, 'weights', weights, overwrite = TRUE)
dbGetQuery(out_conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(out_conn, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(out_conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
dbGetQuery(out_conn, "CREATE INDEX gene_model_summary ON extra (gene)")



library(dplyr)
library(RSQLite)
"%&%" <- function(a,b) paste(a,b, sep='')
tiss <- "AA_hepatocyte"
driver <- dbDriver('SQLite')
# Extra table ----
model_summaries <- read.table('../summary/nested_cv_chr1_model_summaries.txt', header = T, stringsAsFactors = F, fill = T)
tiss_summary <- read.table('../summary/nested_cv_chr1_tiss_chr_summary.txt', header = T, stringsAsFactors = F, fill = T)


for (i in 2:22) {
  model_summaries <- rbind(model_summaries,
                           read.table('../summary/nested_cv_chr' %&% as.character(i) %&% '_model_summaries.txt', header = T, stringsAsFactors = F, fill = T))
  tiss_summary <- rbind(tiss_summary,
                        read.table('../summary/nested_cv_chr' %&% as.character(i) %&% '_tiss_chr_summary.txt', header = T, stringsAsFactors = F, fill = T))
}
#print(model_summaries[1:3,1:3])

model_summaries <- dplyr::rename(model_summaries, gene = gene_id)
write.table(model_summaries, "../dbs/model_summaries.txt", quote =  F, row.names = T, col.names = T)
write.table(tiss_summary, "../dbs/tissue_summaries.txt", quote =  F, row.names = T, col.names = T)
print(dim(model_summaries))

print(model_summaries[1:3,1:3])
conn <- dbConnect(drv = driver, '../dbs/AA_hepatocyte.db')
dbWriteTable(conn, 'model_summaries', model_summaries, overwrite = TRUE)
dbGetQuery(conn, "CREATE INDEX gene_model_summary ON model_summaries (gene)")

# Weights Table -----
weights <- read.table('../weights/nested_cv_chr1_weights.txt', header = T, stringsAsFactors = F, fill = T)
for (i in 2:22) {
  weights <- rbind(weights,
                   read.table('../weights/nested_cv_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F, fill = T))
}
weights <- dplyr::rename(weights, gene = gene_id)
print(dim(weights))
dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

# Sample_info Table ----
sample_info <- data.frame(n_samples = 60, population = 'AA', tissue = "Hepatocyte")
dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)

# Construction Table ----
tiss_summary <- as_data_frame(tiss_summary)
construction <- tiss_summary %>%
  dplyr::select(chrom, cv_seed) %>%
  dplyr::rename(chromosome = chrom)
dbWriteTable(conn, 'construction', construction, overwrite = TRUE)


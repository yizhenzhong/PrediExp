gene <- read.table("../../PrediXcan/data/intermediate/annotations/gene_annotation/gencode.v25lift37.annotation.txt",
                   header = T, stringsAsFactors = F)

gene$id = sapply(gene$gene_id, function(x) strsplit(x, "_")[[1]][1])
exp <- read.table("../../PrediXcan/data/input/expression_phenotypes/baseline_hg19_expr.bed", header = T, 
                  stringsAsFactors = F, check.names = F)

exp$ID= gene$gene_id[match(exp$ID, gene$id)]

exp_t = t(exp[,-c(1:4)])
colnames(exp_t) = exp$ID
exp_t = exp_t[,-which(is.na(colnames(exp_t)))]
write.table(exp_t, "../../PrediXcan/data/input/expression_phenotypes/baseline_hg19_expr.txt",
            quote=F, row.names = T, col.names=NA)


##########################################
cov <- read.table("../../data/baseline_covariant.txt", header = T, stringsAsFactors = F,
                  check.names = F)

cov_t = t(cov[,-1])
colnames(cov_t) = cov$id
write.table(cov_t, "../../PrediXcan/data/input/expression_phenotypes/baseline_hg19_cov.txt",
            quote=F, row.names = T, col.names=NA)



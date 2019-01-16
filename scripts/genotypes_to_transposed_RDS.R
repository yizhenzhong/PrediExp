# Reads an expression file in as a dataframe, transposes it, and
# saves it as an RDS object.  If a covariate file is present (as 3rd
# input argument), then it will be used to produce new expression data
# by making a linear model expression ~ covariate, and then pulling the
# residuals as the new expression data.
#
# The input expression file is expected to be tab-delimted, with people
# as columns and genes as row.

argv <- commandArgs(trailingOnly = TRUE)
genotypefile <- argv[1]
RDSout <- argv[2]
# Presence of covariate file suggests to correct for PEER factors, etc.
#covariatefile <- ifelse(length(argv) == 3, argv[3], NA)

genotype <- read.table(genotypefile, stringsAsFactors = FALSE,
    header = TRUE, row.names = 1, check.names = FALSE)
# Transpose expression.
genotype <- t(genotype)

write.table(genotype, paste0(RDSout, "_geno.txt"), quote=F, row.names = T, col.names=NA)
#saveRDS(expression, RDSout)

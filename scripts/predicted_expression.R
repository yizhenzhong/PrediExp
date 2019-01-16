setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/Yilin/PrediXcan/scripts")
aa <- read.table("../VTE_dosage/VTE_RPKM_AA_predicted_expression.txt", header = T, stringsAsFactors = F)

eur <- read.table("../VTE_dosage/VTE_euro_liver_predicted_expression.txt", header = T,
                  stringsAsFactors = F)

gene_list <- read.table("../../gencode.v25lift37.annotation.txt", header = T, stringsAsFactors = F)
gene_list$gene = sapply(gene_list$gene_id, function(x) strsplit(x, ".", fixed = T)[[1]][1])

gene_aa <- sapply(colnames(aa)[3:ncol(aa)], function(x) strsplit(x, ".", fixed = T)[[1]][1])
gene_eur <- sapply(colnames(eur)[3:ncol(eur)], function(x) strsplit(x, ".", fixed = T)[[1]][1])

aa_t = t(aa[,-c(1:2)])[1:3,1:3]

#################Correlation of predicted expression
common_gene <- intersect(gene_aa, gene_eur)

aa_common <- aa[,-c(1:2)][,match(common_gene, gene_aa)]
eur_common <- eur[,-c(1:2)][,match(common_gene, gene_eur)]


cor = c()
for(i in 1:ncol(aa_common)){
  temp = try(cor(aa_common[,i], eur_common[,i]))
  cor = c(cor, temp)
}

pdf("../figure/correlation_VTE_prediciton.pdf")
hist(cor, main = "correlation of predicted expression of AA and EUR model",
     xlab="correlation")
dev.off()

set2 = RColorBrewer::brewer.pal(8, "Set2")


tiff("../../figure/predicted_expression_compare.tiff", height = 15, width = 15, units = "cm", compression = "lzw", res = 1000)

cor =  as.data.frame(cor)
ggplot(cor) + 
  geom_histogram(aes(x=cor),binwidth = 0.1, fill=set2[8], color="black")+
  xlab("correlation of predicted expression")+
  theme(axis.text=element_text(size=20, colour = "black"),
        axis.title=element_text(size=20, colour = "black", 
                                margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.y = element_text(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_blank(),
        legend.key.size = unit(1.5, 'lines'),
        legend.text=element_text(size=20),
        legend.position="bottom", 
        legend.box = "horizontal",
        legend.key=element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))

dev.off()

length(which(cor>=0.8))
length(which(cor<0.8))

####################
#example 
index = which(cor==min(cor, na.rm = T))
example = cbind(aa_common[,index], eur_common[,index])
colnames(example) = c("AA", "GTEx")
example = as.data.frame(example)
tiff(paste0("../../figure/", common_gene[index],".tiff"), height = 15, 
     width = 15, units = "cm", compression = "lzw", res = 1000)
ggplot(example,aes(y=AA,x=GTEx)) + 
  geom_point(size=2.5,  color="black",pch=21)+
  #geom_smooth(method='lm',lwd=1, color=colors[1])+
  xlab("expression in AAs")+
  ylab("expression in GTExs")+
  geom_abline(slope = 1, intercept = 0, color=colors2[2],lwd=1.3, linetype = 2)+
  theme(axis.text=element_text(size=20, colour = "black"),
        axis.title=element_text(size=20, colour = "black", 
                                margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.y = element_text(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title=element_blank(),
        legend.key.size = unit(1.5, 'lines'),
        legend.text=element_text(size=20),
        legend.position="bottom", 
        legend.box = "horizontal",
        legend.key=element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))

dev.off()
########################
#association test
sample <- read.table("../VTE_dosage/VTE_Discovery.sample", header = F, 
                     stringsAsFactors = F)

all(sample$V1 == aa$FID)
all(sample$V1 == eur$FID)

library("foreach")
library("doParallel")

#table = NULL
cl <- makeCluster(4) # create a cluster with 2 cores
registerDoParallel(cl) # register the cluster
res = foreach(h = 3:ncol(aa),
              .combine = "rbind"
) %dopar% {
  predictor <- as.matrix(cbind(aa[,h], sample[,c(4:ncol(sample)-1)]))
  summary(glm(sample[,ncol(sample)] ~ predictor,
              family = binomial))$coefficients[2,]
}
stopCluster(cl) #

res <- cbind(gene_aa, res)

aa_res <- res
aa_res = as.data.frame(aa_res)
aa_res <- aa_res[order(aa_res[,5]),]
aa_res <- cbind(gene_list$gene_name[match(aa_res[,1], gene_list$gene)], aa_res)
aa_res= cbind(aa_res, p.adjust(aa_res[,6], method = 'fdr'))
colnames(aa_res)[1] = "gene_symbol"
aa_res$cor = NA
aa_res$cor[match(common_gene, aa_res$gene_aa)] = cor$cor
write.table(aa_res, "../VTE_dosage/aa_association_res_RPKM.txt", quote = F, col.names = T, row.names = F,
            sep = "\t")
##########################

#table = NULL
cl <- makeCluster(4) # create a cluster with 2 cores
registerDoParallel(cl) # register the cluster
res = foreach(h = 3:ncol(eur),
              .combine = "rbind"
) %dopar% {
  #print(h)
  
  #for(h in 3:(ncol(aa)-2){
  #  print(h)
  predictor <- as.matrix(cbind(eur[,h], sample[,c(4:ncol(sample)-1)]))
  summary(glm(sample[,ncol(sample)] ~ predictor,
              family = binomial))$coefficients[2,]
}

stopCluster(cl) #

res <- cbind(gene_eur, res)

res_eur <- res
res_eur <- cbind(gene_list$gene_name[match(res_eur[,1], gene_list$gene)], res_eur)
res_eur <- res_eur[order(res_eur[,6]),]
res_eur= cbind(res_eur, p.adjust(res_eur[,6], method = 'fdr'))

write.table(res_eur, "../VTE_dosage/EUR_association_res.txt", quote = F, col.names = T, row.names = F,
            sep = "\t")


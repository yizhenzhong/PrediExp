setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/Hepatocyte_project/Yilin/PrediXcan/scripts")
aa <- read.table("../VTE_dosage/VTE_predicted_expression.txt", header = T, stringsAsFactors = F)

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

length(which(cor>=0.8))
length(which(cor<0.8))

########################
#association test
sample <- read.table("../VTE_dosage/VTE_Discovery.sample", header = F, 
                     stringsAsFactors = F)

all(sample$V1 == aa$FID)
predictor <- as.matrix(cbind(aa[,3], sample[,c(4:ncol(sample)-1)]))

library("foreach")
library("doParallel")

#table = NULL
cl <- makeCluster(4) # create a cluster with 2 cores
registerDoParallel(cl) # register the cluster
res = foreach(h = 3:ncol(aa),
              .combine = "rbind"
) %dopar% {
  #print(h)
  
  #for(h in 3:(ncol(aa)-2){
  #  print(h)
  predictor <- as.matrix(cbind(aa[,h], sample[,c(4:ncol(sample)-1)]))
  summary(glm(sample[,ncol(sample)] ~ predictor,
              family = binomial))$coefficients[2,]
}

stopCluster(cl) #

res <- cbind(gene_aa, res)

aa_res <- res
aa_res <- aa_res[order(aa_res[,5]),]

aa_res <- cbind(gene_list$gene_name[match(aa_res[,1], gene_list$gene)], aa_res)
write.table(aa_res, "../VTE_dosage/aa_association_res.txt", quote = F, col.names = T, row.names = F,
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
res_eur <- res_eur[order(res_eur[,5]),]
write.table(res_eur, "../VTE_dosage/EUR_association_res.txt", quote = F, col.names = T, row.names = F,
            sep = "\t")


########################compare the prediction accuracy
library(dplyr)
library(RSQLite)
driver <- dbDriver("SQLite")
ceu_db <- "../gtex_v7_Liver_imputed_europeans_tw_0.5_signif.db"
ceu_conn <- dbConnect(driver, ceu_db)
aa_db <- "../PredictDB_Pipeline_GTEx_v7/model_training/dbs/AA_hepatocyte_tw0.5_signif.db"
aa_conn <- dbConnect(driver, aa_db)

model_aa <- dbGetQuery(aa_conn, "select * from extra")
write.csv(model_aa, "../AA_hepatocyte_tw0.5_signif.csv",quote = F)

model_eur <- dbGetQuery(ceu_conn, "select * from extra")
write.csv(model_eur, "../gtex_v7_Liver_imputed_europeans_tw_0.5_signif_model_summary.csv",quote = F)

model_aa <- read.csv("../AA_hepatocyte_tw0.5_signif.csv", header = T, stringsAsFactors = F)
model_eur <- read.csv("../gtex_v7_Liver_imputed_europeans_tw_0.5_signif_model_summary.csv", header = T, stringsAsFactors = F)

model_aa$gene_short = sapply(model_aa$gene, function(x) strsplit(x, ".", fixed = T)[[1]][1])
model_eur$gene_short = sapply(model_eur$gene, function(x) strsplit(x, ".", fixed = T)[[1]][1])
merged_model <- merge(model_aa, model_eur, by="gene_short", all = F)

pdf("../figure/prediction_accuracy_rho_ave.pdf")
plot(merged_model$rho_avg.x, merged_model$rho_avg.y, pch=20, xlab = "rho in AA", 
     ylab = "rho in EUR", xlim = c(0.15,0.8), ylim=c(0.15, 0.8))
abline(lm( merged_model$rho_avg.y ~ merged_model$rho_avg.x), col ="orange")
dev.off()


model_summary <- read.table("../PredictDB_Pipeline_GTEx_v7/model_training/dbs/model_summaries.txt",
                            header = T, stringsAsFactors = F)
length(which(model_summary$rho_avg>0.1 & model_summary$zscore_pval <0.05))
model_summary$gene_short = sapply(model_summary$gene, function(x) strsplit(x, ".", fixed = T)[[1]][1])
merged_model <- merge(model_summary, model_eur, by="gene_short", all = F)


index = (which(merged_model$rho_avg.x>0.1 & merged_model$zscore_pval <0.05))
cor.test(merged_model$rho_avg.x[index],merged_model$rho_avg.y[index], method = "spearman")

library(ggplot2)
colors = RColorBrewer::brewer.pal(5, "Accent")
colors3 = RColorBrewer::brewer.pal(6, "Paired")
colors2 = RColorBrewer::brewer.pal(6, "Dark2")

merged_model$group = 0
merged_model$group[index] = "Common"
merged_model$group[-index] = "GTEx only"

tiff("../../figure/predixcan_performance_compare.tiff", height = 15, width = 15, units = "cm", compression = "lzw", res = 1000)

ggplot(merged_model,aes(x=rho_avg.x,y=rho_avg.y,fill = group)) + 
  geom_point(size=2.5,  color="black",pch=21)+
  scale_fill_manual(name = "Type", values = c("GTEx only" = colors3[1], "Common" =colors[5])) +
  #geom_smooth(method='lm',lwd=1, color=colors[1])+
  xlab("rho_average in AA")+
  ylab("rho_average in GTEx EUR")+
  geom_abline(slope = 1, intercept = 0, color=colors2[2],lwd=1.3, linetype = 2)+
  theme(axis.text=element_text(size=20, colour = "black"),
        axis.title=element_text(size=20, colour = "black", 
                                margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.y = element_text(),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
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



###########################

tiff("../../figure/predixcan_n_snps_in_window_compare.tiff", height = 15, width = 15, units = "cm", compression = "lzw", res = 1000)
ggplot(merged_model,aes(x=n_snps_in_window.x,y=n_snps_in_window.y)) + 
  geom_point(size=2.5,  color="black",pch=21)+
  #geom_smooth(method='lm',lwd=1, color=colors[1])+
  xlab("#SNPs in window in AAs")+
  ylab("#SNPs in window in GTExs")+
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


##############################
tiff("../../figure/predixcan_n_snps_in_model_compare.tiff", height = 15, width = 15, units = "cm", compression = "lzw", res = 1000)
ggplot(merged_model,aes(y=n.snps.in.model,x=n_snps_in_model)) + 
  geom_point(size=2.5,  color="black",pch=21)+
  #geom_smooth(method='lm',lwd=1, color=colors[1])+
  xlab("#SNPs in model in AAs")+
  ylab("#SNPs in model in GTExs")+
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

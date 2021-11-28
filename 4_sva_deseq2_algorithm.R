####### SVA DEseq
## we will compare against a null model to estimate any unseen variation in the data
## and quantify as surrogate variables, and then create a new DESeq object with this equation
norm_counts_bulk <- counts(dds_covid, normalized = TRUE)
idx_bulk  <- rowMeans(norm_counts_bulk) > 1
norm_counts_bulk <- norm_counts_bulk[idx_bulk, ]
mod_bulk <- model.matrix(~ 0 + cohort + age_category + Sex, colData(dds_covid))
mod0_bulk<- model.matrix(~ 0 + age_category + Sex , colData(dds_covid))
svseq_1 <- svaseq(norm_counts_bulk , mod_bulk, mod0_bulk)
## 4 surrogate vars                       

### GRAPHS OF SURROGATE VARS #####
par(mfrow = c(4, 2), mar = c(3,3,3,1))
for (i in 1:4) {
  stripchart(svseq_1$sv[, i] ~ dds_covid$cohort, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}
par(mfrow = c(2, 1), mar = c(3,3,3,1))


## create surrogate vars 
ddssva_1 <- dds_covid
ddssva_1$SV1 <- svseq_1$sv[,1]
ddssva_1$SV2 <- svseq_1$sv[,2]
ddssva_1$SV3 <- svseq_1$sv[,3]
ddssva_1$SV4 <- svseq_1$sv[,4]
ddssva_1$SV5 <- svseq_1$sv[,5]

design(ddssva_1) <- ~ 0 + cohort + age_category  + SV1 + SV2 + SV3 + SV4 + SV5

############ DESEQ ############
## Run standard DESEQ2 algorithm on our object
ddssva_1 <- DESeq(ddssva_1)
rowData(ddssva_1)
ressva_1 <- DESeq2::results(ddssva_1, contrast = c("cohort", "COVID", "Control"), alpha = 0.5)
summary(ressva_1)
ressva_lfc <- DESeq2::results(ddssva_1, contrast = c("cohort", "COVID", "Control"), lfcThreshold = 1, alpha = 0.1)
summary(ressva_lfc )



## MA Plots
par(mfrow = c(2, 2))
plotMA(ressva_1 , ylim = c(-8,8))
plotMA(ressva_lfc, ylim = c(-8,8))


## ## list of significants
resSig_1<- subset(ressva_1, padj < 0.05)
head(ressva_1[ order(ressva_1$padj), ])
## go ahead and write csv
write.csv(resSig_1, "/Users/beccabell/RNASeq/https:/github.com/rebell90/Transcriptome_COVID_Control_GSE176480/xx_Data_Files/SVA_sigfiicant_p05.csv")

resSig_lfc<- subset(ressva_lfc, padj < 0.1)
head(ressva_lfc[ order(ressva_lfc$log2FoldChange), ])
## write csv
write.csv(resSig_lfc, "/Users/beccabell/RNASeq/https:/github.com/rebell90/Transcriptome_COVID_Control_GSE176480/xx_Data_Files/SVA_sigfiicant_lfc1_p1.csv")


### in order to examine individual genes of significance, we are going to
### examine machine learning techniques from the DaMiRseq package

###### DESEQ2 ######
## create DeSeq Dataset and preform prelimenary QC Analysis:


# In DESEQ2, the geometric mean is calculated for each gene across all samples. This number is then used as a divisor
# for the gene counts in each sample. The median of the ratios in a sample is then calculated to give a  "size-factor" for
# that specific sample.  This normalizes/corrects for any differences in library size and composition biases. 
# Next, dispersion values are estimated for each gene using a shrinkage model.
# In the DESeq test, a negative binomial linear model is used for each gene and then significance is tested through the  Wald test
# Outliers are also removed using Cook's distance.


## CREATE DESeq OBJECT
dds_covid <- DESeqDataSetFromMatrix(countData = GSE176480_counts, ## count matrix
                                    colData = GSE176480_col, ## sample data
                                    design = ~ 0 + cohort + age_category + Sex) ## design matrix (cohort as the main variable, adjusting for age/sex as covariates)

## pre-filter
keep_feature <- rowSums(counts(dds_covid)  > 20 ) >= 8
dds_covid <-  dds_covid[keep_feature, ]


## box plot
boxplot(log10(counts(dds_covid )+1), las=2) 


### mean-var plots 
## in order to perform certain measures like PCA, the variance across 
## different levels of the mean is expected to be close to the same value (homoskedasticity).  
## in RNA seq data, the variance is expected to grow with the mean, so transformations
## are recommend for any additonal EDA analysis (NOTE: for the actual DESeq algorithm, we will use the raw counts )                        

msd <- meanSdPlot(counts(dds_covid), ranks = FALSE)
msd

## log plus 1                    
log.cts.one <- log2(counts(dds_covid) + 1)
log.cts.one_msd <- meanSdPlot(log.cts.one , ranks = FALSE)
log.cts.one_msd

# NOTE: even the log transform make the variance of LOWER counts to be higher, 
# so, there are other recommended transformations to further stabilize                        

## two options are the VST (Variance-Stabilizing Transformation) or the 
## rlog (log-regularization).
## BOth of these transformations will output similar variances to the log-plus-one
## transforms in regards to larger countss, but for genes with lower counts, the values
## are shrunk towards a middle value 
## This makes the data usable for PCA, distance measures, or any other techiques that have 
## the assumptions of homoskedasticity.



## variance stabilizing with counts
vsd_1 <- DESeq2::vst(dds_covid, blind = TRUE)
head(vsd_1, 3)
## box plot
par(mar=c(5,5,4,2)+0.1)
boxplot(assay(vsd_1), xlab="", ylab="Log2 counts per million",las=2, main="VST Normalized Count Distributions")
abline(h=median(assay(vsd_1)), col="blue")

# rlog
rld_1 <- rlog(dds_covid, blind = TRUE)
head(assay(rld_1), 3)
## box plot
par(mar=c(12,5,4,2)+0.1)
boxplot(assay(rld_1), xlab="", ylab="Log2 counts per million",las=2, main="RLD Normalized Count Distributions")
abline(h=median(assay(rld_1)), col="blue")

### Compare the plots of the 3 transformations
dds_covid <- estimateSizeFactors(dds_covid)

df_1 <- bind_rows(
  as_tibble(log2(counts(dds_covid, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd_1)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld_1)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df_1)[1:2] <- c("x", "y")  
lvls <- c("log2(x + 1)", "vst", "rlog")
df_1$transformation <- factor(df_1$transformation, levels=lvls)
compare_transforms <- ggplot(df_1, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
compare_transforms

##### Distance measures -- see how similar samples are to eachother 

## VSD
sampleDists_vsd <- dist(t(assay(vsd_1)))
sampleDists_vsd

sampleDistMatrix_vsd_1 <- as.matrix(sampleDists_vsd)
rownames(sampleDistMatrix_vsd_1) <- paste(vsd_1$cohort)
colnames(sampleDistMatrix_vsd_1) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "BrBG")))(255)

sampleDists_vsd_1_plot <- pheatmap(sampleDistMatrix_vsd_1,
                                   clustering_distance_rows = sampleDists_vsd,
                                   clustering_distance_cols = sampleDists_vsd,
                                   col = colors)

sampleDists_vsd_1_plot 

## RLOG
sampleDists_rld <- dist(t(assay(rld_1)))
sampleDists_rld

sampleDistMatrix_rld_1 <- as.matrix(sampleDists_rld)
rownames(sampleDistMatrix_rld_1) <- paste(rld_1$cohort)
colnames(sampleDistMatrix_rld_1) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "BrBG")))(255)
sampleDists_rld_1_plot <- pheatmap(sampleDistMatrix_rld_1,
                                   clustering_distance_rows = sampleDists_rld,
                                   clustering_distance_cols = sampleDists_rld,
                                   col = colors)
sampleDists_rld_1_plot 


####### PCA (QUICK) ######
plotPCA(vsd_1, intgroup = "cohort")
plotPCA(rld_1, intgroup = "cohort")

#### ****** SEE sub-section 2b for a more in-depth look and exploration with PCA



### RUN STANDARD DESEQ2 algorithm

dds_covid <- DESeq(dds_covid) ## run deseq2 algorithm

# create results object
res <- DESeq2::results(dds_covid, contrast = c("cohort", "Control", "COVID") , alpha = 0.05) 

### below is the summary of the results object, which displays genes that are up/down-regulated, 
### given the parameters for the results function (i.e. above is for all with p < 0.05)
summary(res) ## summary
resultsNames(dds_covid ) ## names of coefficents (for the cohort variable, there is only one
                          # contrast, given our variable of interest is simpoly binary)

table(res$padj < 0.05) ## 

# MA Plot
plotMA(res, ylim = c(-5, 5))

## with LFC --- here we increase the value of log2 fold change threshold from 0 to 1 , meaning  there 
## has to be at least 2 times as much difference between the COVID and Control group in order for
## a gene to be considered significant
## needed 
res_lfc <- DESeq2::results(dds_covid, contrast = c("cohort", "Control", "COVID") , alpha = 0.1, lfcThreshold=1)
summary(res_lfc)
table(res_lfc$padj < 0.1) ## raise p-val to .1, since the LFC is now so much higher
# FALSE = 8183   TRUE = 149

## list of signficant for lfc
resSig <- subset(res_lfc, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
## MA Plot
plotMA(res_lfc, ylim = c(-5, 5))

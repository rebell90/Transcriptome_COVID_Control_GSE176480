
## PCATOOLS

## explore more visuals/stats on PCA Analysis 
## Will continue to update
p <- pca(assay(vsd_1), metadata = colData(dds_covid), removeVar = 0.1)

## visual depicting the explaned variation percentage for each PCA
screeplot(p, axisLabSize = 18, titleLabSize = 22)

## standard PCA plot of first PCAs , with genes labeled on the eigenvectors, and 
## observations labeled with sample names
biplot(p, showLoadings = TRUE,
       labSize = 3, pointSize = 5, sizeLoadingsNames = 5)

pairsplot(p)

plotloadings(p, labSize = 3)
eigencorplot(p,
             metavars = c('age_category','cohort','Sex','vital_status'))
## we see which PCs help differentiate each of the covariates

## internal data
p$rotated[1:5,1:5]
p$loadings[1:5,1:5]

## PCAs to retain
horn <- parallelPCA(assay(vsd_1))
horn$n
#### 4

elbow <- findElbowPoint(p$variance)
elbow 
#### 3

library(ggplot2)

screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))


#### going to look into this and other dimensionality reductions later! 
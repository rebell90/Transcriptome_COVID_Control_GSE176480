## Assessment of the platelet transcriptome in patients with COVID-19 versus controls ##

## GEO SUMMARY https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
# A cohort of eight hospitalized COVID-19 patients (n=8) were recruited from NYU Langone Health between May 11-21, 2020. 
# SARS-CoV-2 infection was confirmed by RT-PCR, in accordance with current standards. 
# All COVID-19 patients and , age-, and sex-matched control donors were recruited under study protocols approved by the NYU Langone Health Institutional Review Board. 
# Each study participant or their Legal Authorized Representative gave written informed consent for study enrollment in accordance with the Declaration of Helsinki. 
# For COVID-19 patients, enrollment criteria included age greater than 18, hospital admission, positive SARS-CoV-2 testing, and informed consent. 
# COVID-19 patients were monitored until discharge or death.

mem.maxNSize(nsize = Inf)
mem.maxVSize(vsize = Inf)

## download data from GEO
GSE176480 <- getGEO("GSE176480")

## Parse phenoData from GEO
GSE176480 <- GSE176480[["GSE176480_series_matrix.txt.gz"]]
GSE176480_col <- pData(GSE176480)
GSE176480_col 

## Data-Cleaning/Analyzing phenoData 
colnames(GSE176480_col) ## check column names

# format column names into a more R-friendly format
colnames(GSE176480_col) <- gsub(":ch1", "", colnames(GSE176480_col) , fixed = TRUE)
GSE176480_col$title <- gsub("platelets, ", "", GSE176480_col$title , fixed = TRUE)

## select variables to potentially bring into the model
vars <- c("title", "geo_accession", "age_category", "cohort", "Sex", "vital_status" )
GSE176480_col <- GSE176480_col[, vars] # subset

# if the vital status in NULL, it can be assumed that this is one of the control samples, 
# so, we will update this accordingly so the data isnt dropped upon analysis, then convert to factor
GSE176480_col$vital_status <- ifelse(GSE176480_col$vital_status == 'NA', 'Control', GSE176480_col$vital_status)

## convert all character-vars to factor
GSE176480_col <- GSE176480_col %>%  ## convert remaining vars to factors 
  mutate_if(is.character,as.factor) 

## populate row names with the title (these rownames must match the columns in the count-data that will be added)
rownames(GSE176480_col) <- GSE176480_col$title


##### COUNT MATRIX #####

## load data 
GSE176480_counts <- fread("/Users/beccabell/RNASeq/https:/github.com/rebell90/Transcriptome_COVID_Control_GSE176480/xx_Data_Files/GSE176480_covidplatelet_featurecounts.txt")
# view data
head(GSE176480_counts , 20)

## this must be in a format where the gene-identifier is the rowname , and the 
## sample identifier is the column name. The only data should be numeric counts in matrix format

## populate the rownames with the values from the "#GENE" column and then remove the #GENE col
GSE176480_counts <- column_to_rownames(GSE176480_counts, "#GENE")

## subsetting : because these matricies in DGE tend to be very sparse , it is good to eliminate genes
## that have 0 counts in all samples, along with those which such low counts that will likely not show any
## significance.  Not only will this free up memory and storage, it will also lead to more accurate results

## Start by eliminating rows with all 0 counts
keep <- rowSums(GSE176480_counts > 0) > 0 ## save a vector of only rows that have at least a count of 1 in 1 or more genes
GSE176480_counts <- GSE176480_counts[keep,] ## subset the count matrix to genes that meet the criteria above
dim(GSE176480_counts) ## DOUBLE CHECK DIMENSIONS
# 30313    18


## Create a list of thse gene names and save to an object (will use later for annotation data)
gene_list <- rownames(GSE176480_counts)

###### ANNOTATION DATA ######

# We will query data from Annotation DBI, to obtain vital gene metadata, and then we will
# filter this through the count matrix, so we will only be analyzing genes that have a 
# thorough annotation
library(AnnotationDbi)

# Load the Ensembl library
library(EnsDb.Hsapiens.v86)

# Check object metadata
EnsDb.Hsapiens.v86

# Explore the fields that can be used as keys
keytypes(EnsDb.Hsapiens.v86)

# filter the data with the ensembl ids with our count matrix gene-list, using gene-symbol as the matching key
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = gene_list ,
                                         columns = c("GENEID", "ENTREZID","GENEBIOTYPE"),
                                         keytype = "SYMBOL")


library(dplyr)
## group data by symbol, so there is only 1 record per symbol
annotations_edb2 <- annotations_edb %>% 
  dplyr::group_by(SYMBOL) %>% 
  dplyr::summarize(GENEID = paste(GENEID, collapse = ", "),
                   GENEBIOTYPE = paste(GENEBIOTYPE, collapse = ", "))

## populate rownames with the symbol
annotations_edb2 <- column_to_rownames(annotations_edb2, "SYMBOL")

### filter out any genes not in annotation data
## subset counts to remove controls
GSE176480_counts <- GSE176480_counts[which(rownames(GSE176480_counts) %in% rownames(annotations_edb2)),]
## filter out mitochondrial data denoted by "-MT" (these usually have low quality readings)
GSE176480_counts <- GSE176480_counts[-(which(str_sub(rownames(GSE176480_counts), 1, 3) == 'MT-')),]
## order cols/rows (must be matching when the Differential Gene Expression Objects are created)
GSE176480_col <- GSE176480_col[ order(row.names(GSE176480_col)), ]
GSE176480_counts <- GSE176480_counts[, order(colnames(GSE176480_counts)) ]

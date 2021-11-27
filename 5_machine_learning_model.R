
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DaMiRseq/inst/doc/DaMiRseq.pdf

#### Machine Learning approach

## 1) normalization (quantify batch effects, covariates, biological and technical variation, along with any other hidden covariates)
## 2) feature selection to reduce the large number of genes to a small set of powerful predictors
## 3) Classification using a mix of weak classifiers ( aka "Stacked" Generalization method)
## Random Forest (RF), Naïve Bayes (NB), 3-Nearest Neighbours (3kNN), Logistic Regression (LR),
## Linear Discriminant Analysis (LDA), Support Vectors Machines (SVM), Neural Networks (NN) and Partial Least Squares (PLS); \

covariate_cohort <- GSE176480_col %>%
  dplyr::rename(class = cohort) %>%
  dplyr::select(age_category, Sex,  class)

covariate_cohort$class <- droplevels.factor(covariate_cohort$class , exclude = if(anyNA(levels(covariate_cohort$class ))) NULL else NA)
## make summarized experiment object
SE_covid_cohort <-DaMiR.makeSE(GSE176480_counts, covariate_cohort)

## filter out genes with high variability and/or low expression and then normalize
data_norm_cohort <- DaMiR.normalization(SE_covid_cohort, minCounts=80, fSample=0.8,
                                        hyper = "no")

## rerun to remove “hypervariants” -- genes that present anomalous read count-- this is done by comparing to the mean.
## the genes are identified by calculating distinct CV on sample sets that belong to each ’class’ and then removed
data_norm_cohort <- DaMiR.normalization(SE_covid_cohort, minCounts=80, fSample=0.8,
                                        hyper = "yes", th.cv=3)

## samples that are lowly correlated with technical replicates can indicate a technical artifact/error that could throw off later analysis

## these samples are filtered out by assessment of the mean absolute correlation of each sample and removes samples lower than 
## a value set in th.corr argument
data_filt_cohort  <- DaMiR.sampleFilt(data_norm_cohort , th.corr=0.6)

## test for the presence of surrogate variables to remove effects of confounding variables -- the algorithms cannot tell the 
## difference between biological confounders verses technical error variables, so this helps us explore any significant, hidden variation
## contained within a potential confounding variable 

sv_cohort <- DaMiR.SV(data_filt_cohort)

## analysis of whether the surrogate variables are sources of unwanted variation or somehow connected to the other covariates ## looks like 1 and 2 have some correlation with class, but 2 also with gender
DaMiR.corrplot(sv_cohort , colData(data_filt_cohort),  sig.level = 0.01)

## batch effects
data_adjust_cohort <-DaMiR.SVadjust(data_filt_cohort, sv_cohort, n.sv=3)

## QC Plots
DaMiR.Allplot(data_filt_cohort, colData(data_filt_cohort))

#After sample filtering and sv adjusting
DaMiR.Allplot(data_adjust_cohort , colData(data_adjust_cohort))

######## MACHINE LEARNING TASKS #########
## per vingette
# 1) Finding a small set of robust features to discriminate classes.
# 2)  Building the most reliable model to predict new samples.

########### task 1 #########
### subset features (aka genes)
# identify principal components (PCs) that correlate with “class” and then iteratively use backwards selection algorithm
set.seed(12345) ## set seed for reproducability
data_clean_cohort <-DaMiR.transpose(assay(data_adjust_cohort))
df_cohort<-colData(data_adjust_cohort)
data_reduced_cohort <- DaMiR.FSelect(data_clean_cohort, df_cohort, th.corr=0.4)
### 40 genes remain

# remove highly correlated features by creating a pair-wise absolute correlation matrix, and when 2 features are correlated above the 
# pre-defined cutoff, the algorithm calculates mean absolute correlation of each feature and, then, removes the feature with the 
# largest mean absolute correlation

data_reduced_cohort <- DaMiR.FReduct(data_reduced_cohort$data)
## 23 Highly correlated features have been discarded for classification. 
## 15 Features remained. 

## MDS Plot
DaMiR.MDSplot(data_reduced_cohort, df_cohort)

## re calc relative importance of genes
par(margin(1,1))
df.importance_cohort <- DaMiR.FSort(data_reduced_cohort, df_cohort)

## subset further to best predictor-genes: 10 PREDICTORS
selected_features_cohort <- DaMiR.FBest(data_reduced_cohort, ranking = df.importance_cohort, n.pred = 10)

## clust Plot
DaMiR.Clustplot(selected_features_cohort$data, df_cohort)

######## TRAIN Classifier:
Classification_res_cohort <- DaMiR.EnsembleLearning(selected_features_cohort$data,
                                                    classes=df_cohort$class, fSample.tr = 0.7,
                                                    fSample.tr.w = 0.7, iter = 30)


set.seed(10101)
nSampl_cl1_cohort <- 2
nSampl_cl2_cohort <-2
# Create balanced Learning and Test sets
idx_test_cl1_cohort <-sample(1:(ncol(data_adjust_cohort)/2), nSampl_cl1_cohort)
idx_test_cl2_cohort<-sample(1:(ncol(data_adjust_cohort)/2), nSampl_cl2_cohort) + ncol(data_adjust_cohort)/2
idx_test_cohort <- c(idx_test_cl1_cohort, idx_test_cl2_cohort)
Test_set_cohort <- data_adjust_cohort[, idx_test_cohort, drop=FALSE]
Learning_set_cohort <- data_adjust_cohort[, -idx_test_cohort, drop=FALSE]

# Training and Test into a 'nfold' Cross Validation
nfold_cohort <- 3
cv_sample_cohort <- c(rep(seq_len(nfold_cohort), each=ncol(Learning_set_cohort)/(2*nfold_cohort)),
                      rep(seq_len(nfold_cohort), each=ncol(Learning_set_cohort)/(2*nfold_cohort)))
# Variables initialization
cv_models_cohort <- list()
cv_predictors_cohort <- list()
res_df_cohort <- data.frame(matrix(nrow = nfold_cohort, ncol = 7))
colnames(res_df_cohort) <- c("Accuracy",
                             "N.predictors",
                             "MCC",
                             "sensitivity",
                             "Specificty",
                             "PPV",
                             "NPV")

### RUN ####
for (cv_fold_cohort in seq_len(nfold_cohort)){
  # Create Training and Validation Sets
  idx_cv_cohort <- which(cv_sample_cohort != cv_fold_cohort)
  TR_set_cohort <- Learning_set_cohort[,idx_cv_cohort, drop=FALSE]
  Val_set_cohort <- Learning_set_cohort[,-idx_cv_cohort, drop=FALSE]
  #### Feature selection
  data_clean_cohort <-DaMiR.transpose(assay(TR_set_cohort))
  df_cohort <-colData(TR_set_cohort)
  data_reduced_cohort <- DaMiR.FSelect(data_clean_cohort, df_cohort, th.corr=0.6)
  data_reduced_cohort <- DaMiR.FReduct(data_reduced_cohort$data,th.corr = 0.8)
  df_importance_cohort <- DaMiR.FSort(data_reduced_cohort,
                                      as.data.frame(colData(TR_set_cohort)))
  selected_features_cohort <- DaMiR.FBest(data_reduced_cohort,
                                          ranking=df_importance_cohort,
                                          autoselect = "yes")
  TR_set_cohort <- TR_set_cohort[selected_features_cohort$predictors, drop=FALSE]
  Val_set_cohort <- Val_set_cohort[selected_features_cohort$predictors,drop=FALSE]
  ### Model building
  ensl_model_cohort <- DaMiR.EnsL_Train(TR_set_cohort, cl_type = c("RF", "LR"))
  cv_models_cohort[[cv_fold_cohort]] <- ensl_model_cohort
  ### Model testing
  res_Val_cohort <- DaMiR.EnsL_Test(Val_set_cohort,
                                    EnsL_model = ensl_model_cohort)
  # Store all results 
  res_df_cohort[cv_fold_cohort,1] <- res_Val_cohort$accuracy[1] 
  res_df_cohort[cv_fold_cohort,2] <- length(res_Val_cohort$predictors)[1] 
  res_df_cohort[cv_fold_cohort,3] <- res_Val_cohort$MCC[1]
  res_df_cohort[cv_fold_cohort,4] <- res_Val_cohort$sensitivity[1]
  res_df_cohort[cv_fold_cohort,5] <- res_Val_cohort$Specificty[1]
  res_df_cohort[cv_fold_cohort,6] <- res_Val_cohort$PPV[1]
  res_df_cohort[cv_fold_cohort,7] <- res_Val_cohort$NPV[1]
  cv_predictors_cohort[[cv_fold_cohort]] <- res_Val_cohort$predictors
}

##  "E2F1"  "IFI27"
## results
res_df_cohort[,1:5]


idx_best_model_cohort <- DaMiR.ModelSelect(res_df_cohort,
                                           type.sel = "mode",
                                           npred.sel = "min")

res_predict_cohort <- DaMiR.EnsL_Predict(Test_set_cohort,
                                         bestModel = cv_models_cohort[[idx_best_model_cohort]])

cv_predictors_cohort[[idx_best_model_cohort]]

id_classifier_cohort <- 1 
table(colData(Test_set_cohort)$class, res_predict_cohort[,id_classifier_cohort])

# Prediction assessment for Logistic regression 
id_classifier_cohort <- 3 
# Logistic regression 
table(colData(Test_set_cohort)$class, res_predict_cohort[,id_classifier_cohort])


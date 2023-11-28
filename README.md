# scIAE </br> 
## 1. Introduction  
  scIAE is an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data. It can be used to perform feature extraction, identify cell type and predict disease status.
  
  Given gene expression matrix and label (cell type annotation or disease status) of training set and the gene expression matrix of testing set, scIAE can provide the predicted label of testing set. If true label of testing set is given, the evaluation criteria (ACC, MeanF1, and MedF1) can be calculated to evaluate the classification effectiveness of scIAE. If the number of base classifiers is set to one, the dimensionality reduction result of testing set can also be obtained.
 
  scIAE corresponds to the following paper:
  
  Yin, Q., Wang, Y., Guan, J., Ji, G.. scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data, Brief Bioinform, 2022, 23(1): bbab508. https://doi.org/10.1093/bib/bbab508
  
## 2. Installation
Depends: 

       R (>= 4.0.2)   
       Python (>= 3.8.3)

Requirements: 

      keras (>= 2.4.3)
      tensorflow (>=2.3.1)
      library("keras")
      library("parallel")
      library("caret")
      library("e1071")
      library("kknn")
      library("rpart") 
      library("rBayesianOptimization")
  
## 3. Quick start

Run `main.R`. The parameters can be changed as below.

### 3.1 Prepare data

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.5168428. If users want to use their own datasets, the order of cells in gene expression matrix should correspond to that in labels. Rows refer to cells, and columns refer to genes.

      train_data <- as.matrix(read.csv("pancreas_smartseq_data.csv", row.names = 1)) #gene expression matrix of training set (matrix or data.frame, not null)
      train_info <- read.csv("pancreas_smartseq_label.csv", row.names = 1)[, 1] #label of training set (character or integer, not null)
      test_data <- as.matrix(read.csv("pancreas_celseq_data.csv", row.names = 1)) #gene expression matrix of testing set (matrix or data.frame, not null)
      test_info <- read.csv("pancreas_celseq_label.csv", row.names = 1)[, 1] #label of testing set (character or integer, should be provided when calculating ACC, MeanF1, and MedF1)

      
      > train_data[1:5,1:5]
                 GCG       PPY     REG1A    INS       SST
      AZ_A10  2.691009 13.919472 2.770087 5.973918 15.554698
      AZ_A11 16.082961  2.459152 3.445946 2.201366  4.281374
      AZ_A12  5.319866  4.941085 3.733110 1.017788 15.036942
      AZ_A2   4.435744 15.090921 3.333762 2.358494  3.447319
      AZ_A5   8.230754  6.415130 6.945927 2.383204  6.895719
     
     > head(train_info)
      [1] "delta"  "alpha"  "delta"  "gamma"  "ductal" "alpha" 
      
      > test_data[1:5,1:5]
                  REG1A  INS    GCG     CHGB    TM4SF4
      D28.1_1  5.314580   0 10.642062 8.069391 6.143180
      D28.1_13 7.106880   0  5.975432 0.000000 0.000000
      D28.1_15 5.314580   0 10.642062 8.207436 6.979185
      D28.1_17 4.518002   0 10.642062 5.725955 5.657918
      D28.1_2  4.216514   0  6.061108 1.588734 1.588734
      
      > head(test_info)
      [1] "alpha"       "ductal"      "alpha"       "alpha"       "endothelial" "endothelial"

### 3.2 Get overlapping genes (optional)

`get_intersection()` can get overlapping genes between training set and testing set. In this case, the gene expression matrices of training set and testing set should have gene names.
  
      > dim(train_data)
      [1]  2166 10698
      > dim(test_data)
      [1] 2122 6878
      > data_intersection <-get_intersection(train_data,test_data)
      > train_data <- data_intersection[[1]]
      > test_data <- data_intersection[[2]]
      > dim(train_data)
      [1] 2166 4943
      > dim(test_data)
      [1] 2122 4943
 
Note that the data used here is the one from the Hemberg lab, which is different from that we uploaded to Zenodo. The datasets we uploaded to Zenodo were pre-processed, including extracting overlapping genes between training set and testing set.

### 3.3 Cross validation (optional)
`cross_validation()` can perform cross validation for tuning parameters of scIAE, including the number of base classifiers, denoising rate, lambda (regularization parameter), activation functions of hidden layer and output layer, and the encoded dimensions in each stack. Moreover, the function can be used to tune the hyperparameters of base classifiers, including the cost and gamma for SVM, the split criterion for DT, the number of neighbors for kNN, and the number of components for PLSDA. The inputs of `cross_validation()` contain intervals of parameters given by users, training data and corresponding label, and the number of folds for cross validation (default: 5). Then, the function can perform cross validation and return ACC, MeanF1, and MedF1 for each parameter combination. Users can choose the parameters to be used based on their preferences.
  
      > cv_result <- cross_validation(train_data, 
                                      train_info,
                                      t_interval = c(5,10,15), 
                                      denoising_rate_interval = c(0.1,0.2,0.3), 
                                      lambda_interval = c(1e-4,1e-5),
                                      base_classifier = 'SVM',
                                      cost_interval = c(8,16),
                                      gamma_interval = c(1/500,1/1000))
      
### 3.4 Run scIAE
`scIAE()` returns predicted results of testing data. Its inputs are listed below.

      train_data: gene expression matrix of training set (matrix or data.frame, not null)
      train_info: label of training set (character or integer, not null)
      test_data: gene expression matrix of testing set (matrix or data.frame, not null)
      t: number of base classifiers (integer, default: 10)
      denoising_rate: denoising rate in the input layer (numeric, default: 0.2)
      lambda: L1 regularization parameter (numeric, default: 1e-5)
      activation_hidden: activation function used in the hidden layer of each stack (in c('linear','sigmoid','tanh','relu','exponential','softmax'), default: 'sigmoid')
      activation_output: activation function used in the output layer of each stack (in c('linear','sigmoid','tanh','relu','exponential','softmax'), default: 'sigmoid')
      batch_size: batch size in training autoencoder (integer, default: 256)
      learning_rate: learning rate in training autoencoder (numeric, default: 0.001)
      epochs: epochs in training autoencoder (integer, default: 40)
      encoded_1: encoded dimension of stack 1 (integer, default: 1024)
      encoded_2: encoded dimension of stack 2 (integer, default: 128)
      base_classifier: base classifier algorithm (in c('SVM','DT','kNN','PLSDA'), default: 'SVM')
      verbose: if current ensemble is printed or not (logical, default: TRUE)
      cost: cost of constraints violation if base_classifier is 'SVM' (numeric, Default:16)
      gamma: parameter for radial basis if base_classifier is 'SVM' (numeric, Default:1/1000)
      split: split rule if base_classifier is 'DT' (in c('gini','information'), Default: 'information')
      kNN_k: number of neighbors if base_classifier is 'kNN' (integer, Default:5)
      n_components: number of components if base_classifier is 'PLSDA' (integer, Default:10)
      unassigned: if the classifier gives 'unassigned' label or not (logical, Default: FALSE)
      unassigned_threshold: the probability threshold of giving 'unassigned' label (numeric, Default: NA)
      DR_output: if the dimensionality reduction result of testing set is returned or not (logical, Default: TRUE)
 
 Run `scIAE()`, then predicted results can be obtained. 
       
      > scIAE_output <- scIAE (train_data,train_info,test_data)
      > pred_labels <- scIAE_output[['pred_labels']] 
      > head(pred_labels)
      [1] "alpha"       "ductal"      "alpha"       "alpha"       "endothelial" "endothelial"   
      
 If `t=1`, then the dimensionality reduction result of testing set can also be obtained.
 
      > scIAE_output <- scIAE (train_data,train_info,test_data,t=1)
      > DR_result <- scIAE_output[['DR_result']]
      > dim(DR_result)
      [1] 2122  128
 

### 3.5 Evaluate the classification effectiveness
If `test_info` is provided, `evaluate()` calculates the evaluation criteria (ACC, MeanF1, and MedF1).

      > true_labels <- test_info
      > result <- evaluate(true_labels, pred_labels)
      > print(result)
      [1] 0.9161169 0.8339408 0.9531429

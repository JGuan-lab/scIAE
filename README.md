# scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data </br> 
## 1. Introduction  
  scIAE is an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data to identify cell type and predict disease status.
  
  Given gene expression matrix and label (cell type annotation or disease status) of training set and the gene expression matrix of testing set, scIAE can provide the predicted label of testing set. If true label of testing set is given, the evaluation criteria (Acc, MeanF1, and MedF1) can be calculated to evaluate the classification effectiveness of scIAE.
 
  scIAE corresponds to the following paper:
  
  Yin, Q., Wang, Y., Guan, J., Ji, G.. scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data.
  
## 2. Installation
Depends: 

       R (>= 4.0.4)   
       Python (>= 3.8.3)
       keras (>= 2.4.3)
       tensorflow (>=2.3.1)
Requirements: 
      
      library("keras")
      library("clue")
      library("parallel")
      library("caret")
      library("e1071")
      library("kknn")
      library("rpart")     
  
## 3. Quick start

Run `main.R`. The parameters can be changed as below.

### 3.1 Prepare data

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.5168428. If users want to use their own datasets, the order of cells in gene expression matrix should correspond to that in labels. Rows refer to cells, and columns refer to genes.
  
      train_data <- read.csv("pancreas_smartseq_data.csv") #gene expression matrix of training set (matrix or data.frame, not null)
      train_info <- read.csv("pancreas_smartseq_label.csv") #label of training set (character or integer, not null)
      test_data <- read.csv("pancreas_celseq_data.csv") #gene expression matrix of testing set (matrix or data.frame, not null)
      test_info <- read.csv("pancreas_celseq_label.csv") #label of testing set (character or integer, should be provided when calculating Acc, MeanF1, and MedF1)
      
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

### 3.2 Get overlapping genes (Optional)

`get_intersection()` can get overlapping genes of training set and testing set. In that case, the gene expression matrix of training set and testing set should have gene names.
  
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
 
Note that the data used here is the one from Hemberg lab, which is different from what we uploaded to Zenodo. The datasets we uploaded to Zenodo were pre-processed, including extracting overlapping genes of training set and testing set.
      
### 3.3 Run scIAE
`scIAE()` returns predicted results of testing data. Its inputs include:

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
      unassigned: if the classifier gives 'unassigned' label or not, i.e. provides prediction rejection option or not (logical, default: FALSE)
      unassigned_threshold: the probability threshold of giving 'unassigned' label (numeric, default: NA)
      verbose: if current ensemble is printed or not (logical, default: TRUE)
 
 Run `scIAE()`, then predicted results can be obtained. 
        
      > pred_labels <- scIAE (train_data, train_info, test_data) 
      > head(pred_labels)
      [1] "alpha"       "ductal"      "alpha"       "alpha"       "endothelial" "endothelial"   

### 3.4 Evaluate the classification effectiveness
If `test_info` is provided, `evaluate()` calculates the evaluation criteria (Acc, MeanF1, and MedF1).

      > true_labels <- test_info
      > result <- evaluate(true_labels, pred_labels)
      > print(result)
      [1] 0.9161169 0.8339408 0.9531429

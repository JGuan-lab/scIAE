# scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data </br> 
## 1. Introduction  
  scIAE is an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data to identify cell type and predict disease status.
  
  Given gene expression matrix and label (cell type annotation or disease status) of training set and the gene expression matrix of testing set, scIAE can provide the predicted label of testing set. If true label of testing set is given, the evaluation criteria (Acc, MeanF1, MedF1) can be calculated to evaluate the classification effectiveness of scIAE.
 
  scIAE correspond to the following paper:
  
  Yin, Q., Wang, Y., Guan, J., Ji, G.. scIAE: an integrative autoencoder-based ensemble classification framework for single-cell RNA-seq data.
  
## 2. Installation
Depends: 

       R (>= 4.0.4)    
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

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.5168428
  
      train_data <- read.csv("kidney_droplet_data.csv") #gene expression matrix of training set (matrix or data.frame, not null)
      train_info <- read.csv("kidney_droplet_label.csv") #label of training set (character or integer, not null)
      test_data <- read.csv("kidney_facs_data.csv") #gene expression matrix of testing set (matrix or data.frame, not null)
      test_info <- read.csv("kidney_facs_label.csv") #label of training set (character or integer, null)
      
      > train_data[1:5,1:5]  
                                 spink3  slc34a1     miox      kap    aldob
      10X_P4_5_AAACCTGAGATGCCAG 0.000000 1.885036 0.000000 2.876936 1.333235
      10X_P4_5_AAACCTGAGTGTCCAT 1.483756 2.056551 2.056551 2.683615 0.000000
      10X_P4_5_AAACCTGCAAGGCTCC 1.240926 0.000000 0.000000 2.757024 0.000000
      10X_P4_5_AAACCTGTCCTTGCCA 0.000000 0.000000 1.948658 3.370391 0.000000
      10X_P4_5_AAACGGGAGCTGAACG 0.000000 1.463778 2.033989 2.659493 0.000000 
      > head(train_info)  
      [1] "leukocyte"                              "endothelial cell"                       "mesangial cell"                        
      [4] "kidney cell"                            "kidney collecting duct epithelial cell" "epithelial cell of proximal tubule"    
      > test_data[1:5,1:5]
                                spink3   slc34a1  miox    kap   aldob
      A1.B001717.3_38_F.1.1    1.132377 3.22624 3.267144   0  3.7626061
      A10.B002775.3_39_F.1.1   0.000000 0.00000 0.000000   0  0.0000000
      A10.MAA000752.3_10_M.1.1 0.000000 0.00000 0.000000   0  0.0000000
      A11.MAA000801.3_11_M.1.1 0.000000 0.00000 0.000000   0  0.0000000
      A12.B001717.3_38_F.1.1   0.000000 0.00000 0.000000   0  0.6571724
      > head(test_info)
      [1] "macrophage"                             "endothelial cell"                       "endothelial cell"                      
      [4] "kidney collecting duct epithelial cell" "macrophage"                             "kidney collecting duct epithelial cell" 

### 3.2 Run scIAE
`scIAE()` returns predicted results of testing data. Its inputs include:

      train_data: gene expression matrix of training set (matrix or data.frame, not null)
      train_info: label of training set (character or integer, not null)
      test_data:gene expression matrix of testing set (matrix or data.frame, not null)
      t: number of base classifiers (integer, Default: 10)
      denoising_rate: denoising rate in the input layer (numeric, Default: 0.2)
      lambda: L1 regularization rate (numeric, Default:1e-5)
      activation_hidden: activation function used in the hidden layer of each stack (in c('linear','sigmoid','tanh','relu','exponential','softmax'),Default: 'sigmoid')
      activation_output: activation function used in the output layer of each stack (in c('linear','sigmoid','tanh','relu','exponential','softmax'),Default: 'sigmoid')
      batch_size: batch size in training autoencoder (integer, Default: 256)
      learning_rate :learning rate in training autoencoder (numeric, Default : 0.001)
      epochs: epochs in training autoencoder (integer, Default: 40)
      encoded_1: encoded dimension of stack 1 (integer, Default: 1024)
      encoded_2: encoded dimension of stack 2 (integer, Default: 128)
      base_classifier: base classifier algorithm(in c('SVM','DT','kNN','PLSDA'), Default:'SVM')
      unassigned: if the classifier gives 'unassigned' label or not (logical, Default: FALSE)
      unassigned_threshold: the probability threshold of giving 'unassigned' label (numeric, Default: NA)
      verbose: if current ensemble is printed or not (logical, Default: TRUE)
 
 Run `scIAE()`, then predicted results can be obtained. Note that  `test_info` is not necessary. If `test_info` is missing, `scIAE()` can still provide the predicted results.
        
      > pred_labels <- scIAE (train_data,train_info, test_data) 
      > head(pred_labels)
      [1] "B"            "Memory CD4 T" "CD14+ Mono"   "Memory CD4 T" "Naive CD4 T"  "CD14+ Mono" 

### 3.3 Evaluate the classification effectiveness
If `test_info` is provided, `evaluate()` calculates the evaluation criteria (Acc, MeanF1, MedF1).

      > true_labels <- test_info
      > result <- evaluate(true_labels,pred_labels)
      > print(result)
      [1] 0.8362235 0.870293 0.9295775

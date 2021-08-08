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

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.5168428. If users want to use their own datasets,  the order of cells of gene expression matrix should correspond to labels. Row refers to cells, and column refers to genes.
  
      train_data <- read.csv("pancreas_human_data.csv") #gene expression matrix of training set (matrix or data.frame, not null)
      train_info <- read.csv("pancreas_human_label.csv") #label of training set (character or integer, not null)
      test_data <- read.csv("pancreas_mouse_data.csv") #gene expression matrix of testing set (matrix or data.frame, not null)
      test_info <- read.csv("pancreas_mouse_label.csv") #label of training set (character or integer, null)
      
      > train_data[1:5,1:5]
                                    A1BG    A1CF    A2M  A2ML1  A4GALT
      human1_lib1.final_cell_0001    0    1.028709   0     0      0
      human1_lib1.final_cell_0002    0    0.000000   0     0      0
      human1_lib1.final_cell_0003    0    0.000000   0     0      0
      human1_lib1.final_cell_0004    0    0.000000   0     0      0
      human1_lib1.final_cell_0005    0    0.000000   0     0      0
      > head(train_info)  
      [1] "acinar" "acinar" "acinar" "acinar" "acinar" "acinar"
      > test_data[1:5,1:5]
                                     X0610007P14Rik X0610009B22Rik X0610009E02Rik  X0610009L18Rik X0610009O20Rik
      mouse1_lib1.final_cell_0001      0.0000000        0.00000      0.0000000              0      0.0000000
      mouse1_lib1.final_cell_0002      0.9138462        0.00000      0.0000000              0      0.0000000
      mouse1_lib1.final_cell_0003      0.0000000        0.00000      0.5906152              0      0.0000000
      mouse1_lib1.final_cell_0004      0.0000000        0.00000      0.0000000              0      0.6576304
      mouse1_lib1.final_cell_0005      0.7252200        0.72522      0.0000000              0      0.0000000
      > head(test_info)
      [1] "beta"    "ductal"  "delta"   "schwann" "delta"   "beta" 

### 3.2 Get intersection genes (Optional)

`get_intersection()` can get intersection genes of training set and testing set. In that case, the gene expression matrix of training set and testing set should have gene names.
  
      > dim(train_data)
      [1]  8569 20125
      > dim(test_data)
      [1]  1886 14878
      > data_intersection <-get_intersection(train_data,test_data)
      > train_data <- data_intersection[[1]]
      > test_data <- data_intersection[[2]]
      > dim(train_data)
      [1]  8569 12473
      > dim(test_data)
      [1]  1886 12473
 
Note that the data used here is original data from Hemberg-lab, which is different from what we uploaded to Zenodo. The datasets we uploaded to Zenodo are pre-processed, including extracting intersection genes of training set and testing set.
      
### 3.3 Run scIAE
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
      [1] "beta"        "ductal"      "delta"       "endothelial" "delta"       "beta"    

### 3.4 Evaluate the classification effectiveness
If `test_info` is provided, `evaluate()` calculates the evaluation criteria (Acc, MeanF1, MedF1).

      > true_labels <- test_info
      > result <- evaluate(true_labels,pred_labels)
      > print(result)
      [1] 0.8419936 0.5534933 0.6829268

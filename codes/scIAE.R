######Help#####
#train_data: gene expression matrix of training set (matrix or data.frame, not null)
#train_info: label of training set (character or integer, not null)
#test_data:gene expression matrix of testing set (matrix or data.frame, not null)
#t: number of base classifiers (integer, Default: 10)
#denoising_rate: denoising rate in the input layer (numeric, Default: 0.2)
#lambda: L1 regularization rate (numeric, Default:1e-5)
#activation_hidden: activation function used in the hidden layer of each stack (in c('linear','sigmoid','tanh','relu','exponential','softmax'),Default: 'sigmoid')
#activation_output: activation function used in the output layer of each stack (in c('linear','sigmoid','tanh','relu','exponential','softmax'),Default: 'sigmoid')
#batch_size: batch size in training autoencoder (integer, Default: 256)
#learning_rate :learning rate in training autoencoder (numeric, Default : 0.001)
#epochs: epochs in training autoencoder (integer, Default: 40)
#encoded_1: encoded dimension of stack 1 (integer, Default: 1024)
#encoded_2: encoded dimension of stack 2 (integer, Default: 128)
#base_classifier: base classifier algorithm(in c('SVM','DT','kNN','PLSDA'), Default:'SVM')
#cost: cost of constraints violation if base_classifier is 'SVM' (numeric, Default:16)
#gamma: parameter for radial basis if base_classifier is 'SVM' (numeric, Default:1/1000)
#split: split rule if base_classifier is 'DT' (in c('gini','information'), Default: 'information')
#kNN_k: number of neighbors if base_classifier is 'kNN' (integer, Default:5)
#n_components: number of components if base_classifier is 'PLSDA' (integer, Default:10)
#unassigned: if the classifier gives 'unassigned' label or not (logical, Default: FALSE)
#unassigned_threshold: the probability threshold of giving 'unassigned' label (numeric, Default: NA)
#DR_output: if dimensionality result of test data is output or not (logical, Default: TRUE)

scIAE <- function(train_data,
                  train_info,
                  test_data,
                  t = 10, 
                  denoising_rate = 0.2, 
                  lambda = 1e-5,
                  activation_hidden = 'sigmoid',
                  activation_output = 'sigmoid',
                  batch_size = 256,
                  learning_rate = 0.001,
                  epochs = 40,
                  encoded_1 = 1024,
                  encoded_2 = 128,
                  base_classifier = 'SVM',
                  cost = 16,
                  gamma = 1/1000,
                  split = 'information',
                  kNN_k = 5,
                  n_components = 10,
                  unassigned = FALSE,
                  unassigned_threshold = NA,
                  DR_output = TRUE
)  {
  source('RP.R')
  source('DR.R')
  source('classify.R')
  library("parallel")
  
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  
  individual_RP = lapply(1:t, function(k) {
    RP <- RP(train_data,test_data)
    return(RP)
  })
  
  print("Random projection is done!")
  
  X <- list()
  for (i in c(1:t)){
    X[[i]]=list(train_data=individual_RP[[i]][[1]],
                test_data=individual_RP[[i]][[2]],
                denoising_rate = 0.2, 
                lambda = 1e-5,
                activation_hidden = 'sigmoid',
                activation_output = 'sigmoid',
                batch_size = 256,
                learning_rate = 0.001,
                epochs = 40,
                encoded_1 = 1024,
                encoded_2 = 128)
  }
  individual_DR = parallel::parLapply(cl,X, DR)

  print("Dimensionality reduction is done!")
  
  Y <- list()
  for (i in c(1:t)){
    Y[[i]] <- list( microtrain_data=individual_DR[[i]][[1]],
                    train_info=train_info,
                    microtest_data=individual_DR[[i]][[2]],
                    base_classifier = 'SVM',
                    cost = 16,
                    gamma = 1/1000,
                    split = 'information',
                    kNN_k = 5,
                    n_components = 10,
                    unassigned = FALSE,
                    unassigned_threshold = NA)
  }
  individual_classifiers <- parallel::parLapply(cl, Y, classify)
  
  test_prediction_matrix <- matrix(unlist(individual_classifiers),length(individual_classifiers[[1]]),length(individual_classifiers))
  if (t == 1){
    test_prediction_result <- test_prediction_matrix
    if (DR_output == TRUE){
      DR_result <- individual_DR[[t]][[2]]
    }
  }else{
    test_prediction_result <- individual_classifiers[[1]]
    for (ii in 1:nrow(test_prediction_matrix)){
      now_row <- test_prediction_matrix[ii, ]
      timesrow <- table(now_row)
      timesrow <- as.data.frame(timesrow)
      maxidx <- timesrow[which(timesrow$Freq == max(timesrow$Freq)), 1]
      maxidx <- as.character(maxidx)
      test_prediction_result[ii] <- maxidx
    }
  }
  pred_labels <- as.character(test_prediction_result)

  print("Classification is done!")
  
  
  if (t == 1 & DR_output == TRUE){
    return(list(pred_labels = pred_labels,
                DR_result = DR_result))
  }
  else{
    return(list(pred_labels = pred_labels))
  }
  
}


######Help#####
#train_data: gene expression matrix of training set (matrix or data.frame, not null)
#train_info: label of training set (character or integer, not null)
#test_data:gene expression matrix of testing set (matrix or data.frame, not null)
#test_info: label of training set (character or integer, not null)
#directory: working directory (character, not null)
#t: number of base classifiers (integer, Default: 10)
#denoising_rate: denoising rate in the input layer (numeric, Default: 0.2)
#lambda: L1 regularization rate (numeric, Default:1e-5)
#activation_hidden: activation function used in the hidden layer of each stack
#(in c('linear','sigmoid','tanh','relu','exponential','softmax'),Default: 'sigmoid')
#activation_output: activation function used in the output layer of each stack
#(in c('linear','sigmoid','tanh','relu','exponential','softmax'),Default: 'sigmoid')
#batch_size: batch size in training autoencoder (integer, Default: 256)
#learning_rate :learning rate in training autoencoder (numeric, Default : 0.001)
#epochs: epochs in training autoencoder (integer, Default: 40)
#encoded_1: encoded dimension of stack 1 (integer, Default: 1024)
#encoded_2: encoded dimension of stack 2 (integer, Default: 128)
#base_classifier: base classifier algorithm(in c('SVM','DT','kNN','PLSDA'), Default:'SVM')
#unassigned: if the classifier gives 'unassigned' label or not (logical, Default: FALSE)
#unassigned_threshold: the probability threshold of giving 'unassigned' label (numeric, Default: NA)
#verbose: if current ensemble is printed or not (logical, Default: TRUE)

scIAE <- function(train_data,
                  train_info,
                  test_data,
                  test_info,
                  directory,
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
                  unassigned = FALSE,
                  unassigned_threshold = NA,
                  verbose = TRUE
                  )                  
{
  library("keras")
  library("clue")
  library("parallel")
  library("caret")
  library("e1071")
  library("kknn")
  library("rpart")
  
  setwd(directory)
  source("evaluate.R")
  true_labels_path <- 'true_labels_scIAE.csv'
  pred_labels_path <- 'pred_labels_scIAE.csv'
  
  set.seed(1)
  datrows <- rownames(train_data)
  rownames(train_data) <- NULL
  colnames(train_data) <- NULL
  
  num_gene <- ncol(train_data)
  if (num_gene > 8192) {
    final_proj_dim <- 8192
  } 
  else{
    final_proj_dim <- 8192/2
    while (final_proj_dim > 0) {
      if (final_proj_dim < num_gene){
        break
      }
      else{
        final_proj_dim <- final_proj_dim/2
      }
    }
  }
  
  for (k in 1:t){
    if (verbose){
      print('current ensemble:')
      print(k)
    }
    microtest_data <- test_data
    microtrain_data <- train_data 
    if (k != 1){
      rm(random_proj_cols, encoder_input, encoder, decoder_input, decoder, ae_input,
         autoencoder, microtest_layer1, microtest_layer2, microtrain_layer1, microtrain_layer2)
    }
    
    random_proj_cols <- sample(num_gene, size = final_proj_dim, replace = FALSE)
    microtrain_layer1 <- microtrain_data[, random_proj_cols]
    microtest_layer1 <- microtest_data[, random_proj_cols]
    
    keras::k_clear_session()
    
    tns = encoder_input = keras::layer_input(shape = final_proj_dim)
    tns = encoder_denoise = keras::layer_dropout(object = encoder_input, rate = denoising_rate)
    
    tns = keras::layer_dense(tns, units = encoded_1, activation = activation_hidden, activity_regularizer = keras::regularizer_l1(lambda))
    tns = keras::layer_activation_leaky_relu(tns)
    encoder = keras::keras_model(inputs = encoder_input, outputs = tns)
    
    tns = decoder_input = keras::layer_input(shape = encoded_1)
    tns = keras::layer_dense(tns, units = final_proj_dim, activation = activation_output, activity_regularizer = keras::regularizer_l1(lambda))
    decoder = keras::keras_model(inputs = decoder_input, outputs = tns)
    
    tns = ae_input = keras::layer_input(final_proj_dim)
    tns = decoder(encoder(tns))
    autoencoder = keras::keras_model(inputs = ae_input, outputs = tns)
    keras::compile(autoencoder, optimizer = keras::optimizer_adam(lr = learning_rate), loss = 'mean_squared_error')
    
    keras::fit(autoencoder, microtrain_layer1, microtrain_layer1, batch_size = batch_size, epochs = epochs, verbose = 0)
    
    microtrain_layer2 <- predict(encoder, microtrain_layer1, batch_size = batch_size)
    microtest_layer2 <- predict(encoder, microtest_layer1, batch_size = batch_size)
    
    tns = encoder_input = keras::layer_input(shape = encoded_1)
    tns = encoder_denoise = keras::layer_dropout(object = encoder_input, rate = denoising_rate)
    tns = keras::layer_dense(tns, units = encoded_2, activation = activation_hidden, activity_regularizer = keras::regularizer_l1(lambda))
    tns = keras::layer_activation_leaky_relu(tns)
    encoder = keras::keras_model(inputs = encoder_input, outputs = tns)
    
    tns = decoder_input = keras::layer_input(shape = encoded_2)
    tns = keras::layer_dense(tns, units = encoded_1, activation = activation_output, activity_regularizer = keras::regularizer_l1(lambda))
    decoder = keras::keras_model(inputs = decoder_input, outputs = tns)
    
    tns = ae_input= keras::layer_input(encoded_1)
    tns = decoder(encoder(tns))
    autoencoder = keras::keras_model(inputs = ae_input, outputs = tns)
    keras::compile(autoencoder, optimizer = keras::optimizer_adam(lr = learning_rate), loss = 'mean_squared_error')
    
    keras::fit(autoencoder, microtrain_layer2, microtrain_layer2, batch_size = batch_size, epochs = epochs, verbose = 0)
    
    microtrain_data <- predict(encoder, microtrain_layer2, batch_size = batch_size)
    microtest_data <- predict(encoder, microtest_layer2, batch_size = batch_size)
    
    if (unassigned == FALSE){
      
      if (base_classifier == 'SVM'){
        model <- svm(as.factor(train_info)~., data = microtrain_data, cost = 16, gamma = 1/1000)
        test_prediction_vector <- predict(model, microtest_data, type = 'class')
      }
      
      if (base_classifier == 'DT'){
        model <- rpart(train_info~., data = as.data.frame(microtrain_data), method = "class", 
                       parms = list(split = "information"))
        test_prediction_vector <- predict(model, as.data.frame(microtest_data), type = 'class')
      }
      
      if (base_classifier == 'kNN'){
        model <- kknn(as.factor(train_info)~., as.data.frame(microtrain_data), as.data.frame(microtest_data), k = 5)
        test_prediction_vector <- as.character(model[["fitted.values"]])
      }
      
      if (base_classifier == 'PLSDA'){
        model <- plsda(microtrain_data, as.factor(train_info), ncomp = 10, type = 'class', probMethod = 'softmax')
        test_prediction_vector <- predict(model, microtest_data, type = 'class')
      }
    }
    else{
      if (base_classifier == 'SVM'){
        model <- svm(as.factor(train_info)~., data = microtrain_data, cost = 16, gamma = 1/1000)
        test_prediction_prob <- predict(model, microtest_data, probability = TRUE)
        pred_prob <- attr(test_prediction_prob, "probabilities")
        test_prediction_vector <- c(1:nrow(pred_prob))
        for (j in 1:nrow(pred_prob)){
          now_prob <- pred_prob[j,]
          maxprob <- max(now_prob)
          if (maxprob > unassigned_threshold){
            test_prediction_vector[j] <- names(which(now_prob == max(now_prob)))
          }
          else
          {
            test_prediction_vector[j] <- 'unassigned'
          }
        }
      }
      
      if (base_classifier == 'DT'){
        model <- rpart(train_info~., data = as.data.frame(microtrain_data), method = 'class', 
                       parms = list(split = "information"))
        test_prediction_prob <- predict(model, as.data.frame(microtest_data), probability = TRUE)
        pred_prob <- test_prediction_prob
        test_prediction_vector <- c(1:nrow(pred_prob))
        for (j in 1:nrow(pred_prob)){
          now_prob <- pred_prob[j,]
          maxprob <- max(now_prob)
          if (maxprob > unassigned_threshold){
            test_prediction_vector[j] <- names(which(now_prob == max(now_prob)))
          }
          else
          {
            test_prediction_vector[j] <- 'unassigned'
          }
        }
      }
      
      if (base_classifier == 'kNN'){
        stop("Base classifier kNN cannot be labelled 'unassigned'.")
      }
      
      if (base_classifier == 'PLSDA'){
        model <- plsda(microtrain_data, as.factor(train_info), ncomp = 10, type = 'class', probMethod = 'softmax')
        test_prediction_prob <- predict(model, as.data.frame(microtest_data), type='prob')
        pred_prob <- test_prediction_prob[,,1]
        test_prediction_vector <- c(1:nrow(pred_prob))
        for (j in 1:nrow(pred_prob)){
          now_prob <- pred_prob[j,]
          maxprob <- max(now_prob)
          if (maxprob > unassigned_threshold){
            test_prediction_vector[j] <- names(which(now_prob == max(now_prob)))
          }
          else
          {
            test_prediction_vector[j] <- 'unassigned'
          }
        }
      }
    }
    
    
    test_prediction_vector <- as.character(test_prediction_vector)
    if (k == 1){
      test_prediction_matrix <- test_prediction_vector
    }
    else{
      test_prediction_matrix <- cbind(test_prediction_matrix, test_prediction_vector)
    }
  }
  

  
  if (t == 1)
  {
    test_prediction_result <- test_prediction_matrix
  }
  else{
    test_prediction_result <- test_info
    for (ii in 1:nrow(test_prediction_matrix)){
      now_row <- test_prediction_matrix[ii, ]
      timesrow <- table(now_row)
      timesrow <- as.data.frame(timesrow)
      maxidx <- timesrow[which(timesrow$Freq == max(timesrow$Freq)), 1]
      maxidx <- as.character(maxidx)
      test_prediction_result[ii]= maxidx
    }
  }

  pred_labels <- as.character(test_prediction_result)
  true_labels <- as.character(test_info)
  
  write.csv(true_labels,'true_labels_scIAE.csv',row.names = FALSE)
  write.csv(pred_labels,'pred_labels_scIAE.csv',row.names = FALSE)

  result <- evaluate(true_labels_path, pred_labels_path)
  evaluate_index <- c(result$Acc, mean(result$F1), median(result$F1))

  return(evaluate_index)
}

classify <- function (para=list(microtrain_data,
                                train_info,
                                microtest_data,
                                base_classifier = 'SVM',
                                cost = 16,
                                gamma = 1/1000,
                                split = 'information',
                                kNN_k = 5,
                                n_components = 10,
                                unassigned = FALSE,
                                unassigned_threshold = NA)  )   {

  library("caret")
  library("e1071")
  library("kknn")
  library("rpart")
  
  microtrain_data <- para[["microtrain_data"]]
  train_info <- para[["train_info"]]
  microtest_data <- para[["microtest_data"]]
  base_classifier <- para[["base_classifier"]]
  cost <- para[["cost"]]
  gamma <- para[["gamma"]]
  split <- para[["split"]]
  kNN_k <- para[["kNN_k"]]
  n_components <- para[["n_components"]]
  unassigned <- para[["unassigned"]]
  unassigned_threshold <- para[["unassigned_threshold"]]
  if (unassigned == FALSE){
    
    if (base_classifier == 'SVM'){
      model <- svm(as.factor(train_info)~., data = microtrain_data, cost = cost, gamma = gamma)
      test_prediction_vector <- predict(model, microtest_data, type = 'class')
    }
    
    if (base_classifier == 'DT'){
      model <- rpart(train_info~., data = as.data.frame(microtrain_data), method = "class", 
                     parms = list(split = split))
      test_prediction_vector <- predict(model, as.data.frame(microtest_data), type = 'class')
    }
    
    if (base_classifier == 'kNN'){
      model <- kknn(as.factor(train_info)~., as.data.frame(microtrain_data), as.data.frame(microtest_data), k = kNN_k)
      test_prediction_vector <- as.character(model[["fitted.values"]])
    }
    
    if (base_classifier == 'PLSDA'){
      model <- plsda(microtrain_data, as.factor(train_info), ncomp = n_components, type = 'class', probMethod = 'softmax')
      test_prediction_vector <- predict(model, microtest_data, type = 'class')
    }
  }
  else{
    if (base_classifier == 'SVM'){
      model <- svm(as.factor(train_info)~., data = microtrain_data, cost = cost, gamma = gamma)
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
                     parms = list(split = split))
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
      stop("cannot provide prediction rejection. Parameter 'unassigned' should be set to TRUE.")
    }
    
    if (base_classifier == 'PLSDA'){
      model <- plsda(microtrain_data, as.factor(train_info), ncomp = n_components, type = 'class', probMethod = 'softmax')
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
  return(test_prediction_vector)
}
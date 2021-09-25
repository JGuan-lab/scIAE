cross_validation <- function(data, 
                             label,
                             n_folds = 5,
                             t_interval = 10, 
                             denoising_rate_interval = 0.2, 
                             lambda_interval = 1e-5,
                             activation_hidden_interval = 'sigmoid', 
                             activation_output_interval = 'sigmoid',
                             encoded_1_interval = 1024, 
                             encoded_2_interval = 128, 
                             base_classifier = 'SVM',
                             cost_interval = 16,
                             gamma_interval = 1/1000,
                             split_interval = 'information',
                             kNN_k_interval = 5,
                             n_components_interval = 10){
  source("scIAE.R")
  source("evaluate.R")
  library(rBayesianOptimization)
  folds <- KFold(label, nfolds = n_folds, stratified = TRUE)
  
  for (i in 1:n_folds){
    if (base_classifier == 'SVM'){
      result_each <- matrix(nrow = 0, ncol=13)
    }
    else{
      result_each <- matrix(nrow = 0, ncol=12)
    }
    
    count=0
    test_data=data[folds[[i]],]
    train_data=data[-folds[[i]],]
    test_info=label[folds[[i]]]
    train_info=label[-folds[[i]]]
    for (t in t_interval){
      for (denoising_rate in denoising_rate_interval){
        for (lambda in lambda_interval){
          for (activation_hidden in activation_hidden_interval){
            for (activation_output in activation_output_interval){
              for (encoded_1 in encoded_1_interval){
                for (encoded_2 in encoded_2_interval){
                  if (base_classifier == 'SVM'){
                    for (cost in cost_interval){
                      for (gamma in gamma_interval){
                        count=count+1
                        pred_labels <- scIAE (train_data,
                                              train_info,
                                              test_data,
                                              t=t,
                                              denoising_rate=denoising_rate,
                                              lambda=lambda,
                                              activation_hidden = activation_hidden,
                                              activation_output = activation_output,
                                              encoded_1 = encoded_1,
                                              encoded_2 = encoded_2,
                                              base_classifier = base_classifier,
                                              cost = cost,
                                              gamma = gamma) 
                        true_labels <- test_info
                        metrics <- evaluate(true_labels,pred_labels)
                        result_each=rbind(result_each,c(t,denoising_rate,lambda,activation_hidden,activation_output,
                                                        encoded_1,encoded_2,base_classifier,cost,gamma,
                                                        metrics[1],metrics[2],metrics[3]))
                      }
                    }
                  }
                  else if (base_classifier == 'DT'){
                    for (split in split_interval){
                      count=count+1
                      pred_labels <- scIAE (train_data,
                                            train_info,
                                            test_data,
                                            t=t,
                                            denoising_rate=denoising_rate,
                                            lambda=lambda,
                                            activation_hidden = activation_hidden,
                                            activation_output = activation_output,
                                            encoded_1 = encoded_1,
                                            encoded_2 = encoded_2,
                                            base_classifier = base_classifier,
                                            split = split) 
                      true_labels <- test_info
                      metrics <- evaluate(true_labels,pred_labels)
                      result_each=rbind(result_each,c(t,denoising_rate,lambda,activation_hidden,activation_output,
                                                      encoded_1,encoded_2,base_classifier,split,
                                                      metrics[1],metrics[2],metrics[3]))
                    }
                  }
                  else if (base_classifier == 'kNN'){
                    for (kNN_k in kNN_k_interval){
                      count=count+1
                      pred_labels <- scIAE (train_data,
                                            train_info,
                                            test_data,
                                            t=t,
                                            denoising_rate=denoising_rate,
                                            lambda=lambda,
                                            activation_hidden = activation_hidden,
                                            activation_output = activation_output,
                                            encoded_1 = encoded_1,
                                            encoded_2 = encoded_2,
                                            base_classifier = base_classifier,
                                            kNN_k = kNN_k) 
                      true_labels <- test_info
                      metrics <- evaluate(true_labels,pred_labels)
                      result_each=rbind(result_each,c(t,denoising_rate,lambda,activation_hidden,activation_output,
                                                      encoded_1,encoded_2,base_classifier,kNN_k,
                                                      metrics[1],metrics[2],metrics[3]))
                    }
                  }
                  else if (base_classifier == 'PLSDA'){
                    for (n_components in n_components_interval){
                      count=count+1
                      pred_labels <- scIAE (train_data,
                                            train_info,
                                            test_data,
                                            t=t,
                                            denoising_rate=denoising_rate,
                                            lambda=lambda,
                                            activation_hidden = activation_hidden,
                                            activation_output = activation_output,
                                            encoded_1 = encoded_1,
                                            encoded_2 = encoded_2,
                                            base_classifier = base_classifier,
                                            n_components = n_components) 
                      true_labels <- test_info
                      metrics <- evaluate(true_labels,pred_labels)
                      result_each=rbind(result_each,c(t,denoising_rate,lambda,activation_hidden,activation_output,
                                                      encoded_1,encoded_2,base_classifier,n_components,
                                                      metrics[1],metrics[2],metrics[3]))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    if (base_classifier == 'SVM'){
      if (i==1){
        result_fold <- apply(result_each[,11:13],2,as.numeric)
      }else{
        result_fold <- result_fold+apply(result_each[,11:13],2,as.numeric)
      }
    }
    else{
      if (i==1){
        result_fold <- apply(result_each[,10:12],2,as.numeric)
      }else{
        result_fold <- result_fold+apply(result_each[,10:12],2,as.numeric)
      }
    }
    
  }
  result_fold <- result_fold/n_folds
  cross_validation_result <- cbind(result_each[,1:(ncol(result_each)-3)],result_fold)    
  if (base_classifier == 'SVM'){
    colnames(cross_validation_result)=c('t','denoising_rate','lambda','activation_hidden','activation_output',
                                        'encoded_1','encoded_2','base_classifier','cost','gamma',
                                        'ACC','MeanF1','MedF1')
  }else if(base_classifier == 'DT'){
    colnames(cross_validation_result)=c('t','denoising_rate','lambda','activation_hidden','activation_output',
                                        'encoded_1','encoded_2','base_classifier','split',
                                        'ACC','MeanF1','MedF1')
  }else if (base_classifier == 'kNN'){
    colnames(cross_validation_result)=c('t','denoising_rate','lambda','activation_hidden','activation_output',
                                        'encoded_1','encoded_2','base_classifier','k_kNN',
                                        'ACC','MeanF1','MedF1')
  }else if(base_classifier == 'PLSDA'){
    colnames(cross_validation_result)=c('t','denoising_rate','lambda','activation_hidden','activation_output',
                                        'encoded_1','encoded_2','base_classifier','n_components',
                                        'ACC','MeanF1','MedF1')
  }
  return(cross_validation_result)
}
  
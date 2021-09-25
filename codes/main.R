source("scIAE.R")
source("evaluate.R")
source('get_intersection.R')
source('cross_validation.R')

######   Prepare data   ######
train_data <- read.csv("XXXX.csv")
train_info <- read.csv("XXXX.csv")
test_data <- read.csv("XXXX.csv")
test_info <- read.csv("XXXX.csv")

######   Get intersection genes (optional)   ######
# data_intersection <-get_intersection(train_data,test_data)
# train_data <- data_intersection[[1]]
# test_data <- data_intersection[[2]]

######   Cross validation (optional)   ######
# cv_result <- cross_validation(train_data, train_info)

######   Run scIAE   ######
scIAE_output <- scIAE (train_data,
                       train_info,
                       test_data)

pred_labels <- scIAE_output[['pred_labels']]
DR_result <- scIAE_output[['DR_result']] 

######   Evaluate   ######
true_labels <- test_info
result <- evaluate(true_labels,pred_labels)
print(result)

cv_result <- cross_validation(train_data, 
                              train_info,
                              t_interval = c(10,15,20), 
                              denoising_rate_interval = c(0.1,0.2,0.3), 
                              lambda_interval = c(1e-4,1e-5),
                              base_classifier = 'DT',
                              split_interval = c('information','gini'))


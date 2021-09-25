source("scIAE.R")
source("evaluate.R")
source('get_intersection.R')

# ######   Prepare data   ######
train_data <- read.csv("XXXX.csv")
train_info <- read.csv("XXXX.csv")
test_data <- read.csv("XXXX.csv")
test_info <- read.csv("XXXX.csv")

#######   Get intersection genes (Optional)   ######
# data_intersection <-get_intersection(train_data,test_data)
# train_data <- data_intersection[[1]]
# test_data <- data_intersection[[2]]



######   Run scIAE   ######
scIAE_output <- scIAE (train_data,
                       train_info,
                       test_data,t=1)

pred_labels <- scIAE_output[['pred_labels']]
DR_result <- scIAE_output[['DR_result']] 

######   Evaluate   ######
true_labels <- test_info
result <- evaluate(true_labels,pred_labels)
print(result)




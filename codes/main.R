source("scIAE.R")
source("evaluate.R")
source('get_intersection.R')

######   Prepare data   ######
train_data <- read.csv("XXXX.csv", row.names = 1)
train_info <- read.csv("XXXX.csv")
test_data <- read.csv("XXXX.csv", row.names = 1)
test_info <- read.csv("XXXX.csv")
train_info <- as.character(train_info$train_label)
test_info <- as.character(test_info$test_label)

######   Get intersection genes (Optional)   ######
#data_intersection <-get_intersection(train_data,test_data)
#train_data <- data_intersection[[1]]
#test_data <- data_intersection[[2]]

######   Run scIAE   ######
pred_labels <- scIAE (train_data,
                      train_info,
                      test_data) 

######   Evaluate   ######
true_labels <- test_info
result <- evaluate(true_labels,pred_labels)
print(result)

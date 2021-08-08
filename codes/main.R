source("scIAE.R")
source("evaluate.R")
train_data <- read.csv("XXXX.csv")
train_info <- read.csv("XXXX.csv")
test_data <- read.csv("XXXX.csv")
test_info <- read.csv("XXXX.csv")
data_intersection <-get_intersection(train_data,test_data)
train_data <- data_intersection[[1]]
test_data <- data_intersection[[2]]
pred_labels <- scIAE (train_data,
                      train_info,
                      test_data) 
true_labels <- test_info
result <- evaluate(true_labels,pred_labels)
print(result)

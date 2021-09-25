source('cross_validation.R')

train_data <- read.csv("XXXX.csv")
train_info <- read.csv("XXXX.csv")

result <- cross_validation(train_data, train_info)

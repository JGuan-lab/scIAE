source("scIAE.R")
load("XXXX.RData")
directory="XXXX"
result=scIAE (train_data,
       train_info,
       test_data,
       test_info,
       directory)     
print(result)

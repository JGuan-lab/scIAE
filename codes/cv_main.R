source('cross_validation.R')
load("G:/Yinqingyang/result/train_test/PBMC/data.RData")
result <- cross_validation(train_data, 
                           train_info,
                           n_folds = 2,
                           t_interval = 1,
                           encoded_1_interval = c(512,1024),
                           base_classifier = 'PLSDA',
                           n_components_interval = c(3,5))

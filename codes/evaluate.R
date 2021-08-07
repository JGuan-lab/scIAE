evaluate <- function(true_labels_path,pred_labels_path){
  
  true_labels <- read.csv(true_labels_path)
  pred_labels <- read.csv(pred_labels_path)
  true_labels <- as.character(unlist(true_labels))
  pred_labels <- as.character(unlist(pred_labels))
  
  true_celltypes <- unique(true_labels)
  true_labels <- true_labels[pred_labels != 'unassigned']
  pred_labels <- pred_labels[pred_labels != 'unassigned']
  pred_celltypes <- unique(pred_labels)
  
  TP <- c()
  FP <- c()
  TN <- c()
  FN <- c()
  P <- c()
  R <- c()
  F1 <- c()
  sum_correct <- 0
  
  for (i in c(1:length(true_celltypes))){
    
    TP[i] <- sum((true_labels == true_celltypes[i]) & (pred_labels == true_celltypes[i]))
    FP[i] <- sum((true_labels != true_celltypes[i]) & (pred_labels == true_celltypes[i]))
    FN[i] <- sum((true_labels == true_celltypes[i]) & (pred_labels != true_celltypes[i]))
    TN[i] <- sum((true_labels != true_celltypes[i]) & (pred_labels != true_celltypes[i]))
    
    
    if (TP[i]+FP[i] == 0){
      P[i] <- 0
    }else{
      P[i] <- TP[i] / (TP[i]+FP[i])
    }
    if ((TP[i]+FN[i]) == 0){
      R[i] <- 0
    }else{
      R[i] <- TP[i] / (TP[i]+FN[i])
    }
    if (P[i] == 0 | R[i] == 0){
      F1[i] <- 0
    }else{
        F1[i] <- 2*P[i]*R[i] / (P[i] + R[i])
    }
    sum_correct <- sum_correct + TP[i]
  }
  
  Acc <- sum_correct / length(true_labels)
  
  names(F1) <- true_celltypes
  
  result <- list(Acc = Acc , F1 = F1)
  
}

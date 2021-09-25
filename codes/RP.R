RP <- function(train_data,test_data){
  datrows <- rownames(train_data)
  rownames(train_data) <- NULL
  colnames(train_data) <- NULL
  
  num_gene <- ncol(train_data)
  if (num_gene > 8192) {
    final_proj_dim <- 8192
  } 
  else{
    final_proj_dim <- 8192/2
    while (final_proj_dim > 0) {
      if (final_proj_dim < num_gene){
        break
      }
      else{
        final_proj_dim <- final_proj_dim/2
      }
    }
  }
  
  random_proj_cols <- sample(num_gene, size = final_proj_dim, replace = FALSE)
  traindata_rp <- train_data[, random_proj_cols]
  testdata_rp <- test_data[, random_proj_cols]
  return(list(traindata_rp,testdata_rp))
}
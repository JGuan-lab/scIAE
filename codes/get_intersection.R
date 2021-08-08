get_intersection <- function(train_data,test_data){
  train_genes <- colnames(train_data)
  test_genes <- colnames(test_data)
  train_genes2 <- tolower(train_genes)
  test_genes2 <- tolower(test_genes)
  colnames(train_data) <- train_genes2
  colnames(test_data) <- test_genes2
  xgenes <- intersect(train_genes2,test_genes2)
  train_data <- train_data[,xgenes]
  test_data <- test_data[,xgenes]
  list_data=list(train_data,test_data)
  return(list_data)
}

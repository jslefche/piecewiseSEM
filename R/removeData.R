#' Remove data from the model list
#' 
#' formulas = 0, keep everything
#' formulas = 1, remove all formulas including correlated errors
#' formulas = 2, remove only formula but keep correlated errors
#' formulas = 3, remove correlated errors but keep formula
#' 
#' @keywords export
#' 
removeData <- function(modelList, formulas = 0) {
  
  remove <- c("character", "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")
  
  if(formulas == 1) remove <- c(remove, "formula", "formula.cerror")
  
  if(formulas == 2) remove <- c(remove, "formula")
  
  if(formulas == 3) remove <- c(remove, "formula.cerror")
  
  modelList[!sapply(modelList, function(x) any(class(x) %in% remove))]
  
}
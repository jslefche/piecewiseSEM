#' Multigroup Analysis for Piecewise SEM
#' 
#' @param modelList
#' @param data
#' @param group
#' 
#' @return 
#' 
multigroup <- function(modelList, data, group) {
  
  # refit models to each group
  modelList.list <- lapply(levels(data[, group]), function(i) update(modelList, data = data[data[, group] == i, ]))
  
  names(modelList.list) <- levels(data[, group])

  # conduct LRT
  
  lapply(2:modelList, function(i) anova(modelList.list[[i]], modelList.list[[i - 1]]) )
  
  
  # fit interactions
  
  
  # get coefficients
  lapply(modelList.list, coefs)
}
  
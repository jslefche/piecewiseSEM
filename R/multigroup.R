#' Multigroup Analysis for Piecewise SEM
#' 
#' @param modelList
#' @param group
#' 
#' @return 
#' 
multigroup <- function(modelList, group) {
  
  data <- modelList$data
  
  modelList <- removeData(modelList, formulas = 1)
  
  # refit model with group-interaction
  newModelList <- lapply(modelList, function(i) {
    
    rhs1 <- paste(paste(all.vars_trans(i)[-1], collapse = " + "))
    
    rhs2 <- paste(paste(all.vars_trans(i)[-1], ":", group), collapse = " + ")

    update(i, formula(paste(". ~ ", paste(rhs1, " + ", rhs2))))
    
  } )
  
  # capture output and assign to each group
  coefList <- lapply(levels(data[, group]), function(i) update(modelList, data = data[data[, group] == i, ]))
  
  
  # test for significant interactions
  
  anovaTable <- anova(as.psem(newModelList))[[1]]
  
  anovaInts <- anovaTable[grepl(":", anovaTable$Predictor), ]
  
  global <- anovaInts[anovaInts$P.Value >= 0.05, c("Predictor", "Response")]
  
  # alter output for each group if interaction is non-significant (Assign global value)
  
  # plus have column for "constrained" vs "free" ??
  
}
  

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
  intModelList <- lapply(modelList, function(i) {
    
    rhs1 <- paste(paste(all.vars_trans(i)[-1], collapse = " + "))
    
    rhs2 <- paste(paste(all.vars_trans(i)[-1], ":", group), collapse = " + ")

    update(i, formula(paste(". ~ ", paste(rhs1, " + ", rhs2))))
    
  } )
  
  # capture output and assign to each group
  coefsList <- lapply(unique(data[, group]), function(i) {
    
    m <- update(as.psem(modelList), data = data[data[, group] == i, ])
    
    coefs(m)
    
  } )
  
  names(coefsList) <- unique(data[, group])
  
  coefTable <- coefs(modelList)
 
  # test for significant interactions
  
  anovaTable <- anova(as.psem(newModelList))[[1]]
  
  anovaInts <- anovaTable[grepl(":", anovaTable$Predictor), ]
  
  global <- anovaInts[anovaInts$P.Value >= 0.05, c("Response", "Predictor")]
  
  global$Predictor <- gsub(paste0("(.*):|", group, "(.*)"), "\\1", global$Predictor)
  
  # alter output for each group if interaction is non-significant (Assign global value)
  newCoefsList <- lapply(coefsList, function(i) {
    
    i[which(i$Response %in% global$Response & i$Predictor %in% global$Predictor), ] <-
      coefTable[which(coefTable$Response %in% global$Response & coefTable$Predictor %in% global$Predictor), ]
  
    i[, ncol(i)] <- ifelse(i$Response %in% global$Response & i$Predictor %in% global$Predictor, "c", "")
    
    i[, ncol(i) + 1] <- isSig(i$P.Value)
    
    names(i)[(ncol(i) - 1):ncol(i)] <- ""
    
    return(i) 
    
  } )

  ret <- list(
    group = group,
    global = global,
    anovaInts = anovaInts,
    group.coefs = newCoefsList
  )
  
  class(ret) <- "multigroup.psem"
  
  return(ret)
  
}

#' Print multigroup
#' 
print.multigroup.psem <- function(x, ...) {
  
  name <- deparse(substitute(object))
  
  cat("\nStructural Equation Model of", name, "\n")
  
  cat("\nGroups =", x$group, ":", paste(names(x$group.coefs), collapse = " "))
  
  cat("\n")
  
  cat("\nModel-wide Interactions:\n")
  
  cat("\n", captureTable(x$anovaInts))
  
  for(i in 1:nrow(x$global)) cat("\n", paste(x$global[i, "Predictor"], "->", x$global[i, "Response"], "constrained to the global model"), "\n")
  
  cat("\n---\n\n")
  
  for(i in names(x$group.coefs)) {
    
    cat(paste("Group:", i), "coefficients \n")
    
    cat("\n", captureTable(x$group.coefs[[i]]), "\n")
    
  }
  
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05    c = constrained")
  
  cat("\n---\n\n")

  invisible(x)
  
}


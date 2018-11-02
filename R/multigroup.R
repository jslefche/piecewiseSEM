#' Multigroup Analysis for Piecewise SEM
#' 
#' @param modelList a list of structural equations
#' @param group the name of the grouping variable in quotes
#' 
#' @author Jon Lefcheck <LefcheckJ@@si.edu>
#' 
#' @examples
#' data(meadows)
#' 
#' jutila <- psem(
#' lm(rich ~ elev + mass, data = meadows),
#' lm(mass ~ elev, data = meadows)
#' )
#' 
#' jutila.multigroup <- multigroup(jutila, group = "grazed")
#' 
#' jutila.multigroup
#' 
#' @export
#' 
multigroup <- function(modelList, group) {
  
  name <- deparse(match.call()$modelList)
  
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
 
  coefTable[, ncol(coefTable)] <- "c"
  
  coefTable[, ncol(coefTable) + 1] <- isSig(coefTable$P.Value)
  
  names(coefTable)[(ncol(coefTable) - 1):ncol(coefTable)] <- ""
  
  anovaTable <- anova(as.psem(intModelList))[[1]]
  
  anovaInts <- anovaTable[grepl(":", anovaTable$Predictor), ]
  
  global <- anovaInts[anovaInts$P.Value >= 0.05, c("Response", "Predictor")]
  
  global$Predictor <- gsub(paste0("(.*):|", group, "(.*)"), "\\1", global$Predictor)
  
  if(nrow(global) == nrow(anovaInts)) newCoefsList <- list(global = coefTable) else {
  
    # alter output for each group if interaction is non-significant (Assign global value)
    newCoefsList <- lapply(coefsList, function(i) {
      
      suppressWarnings(
        i[which(i$Response %in% global$Response & i$Predictor %in% global$Predictor), ] <-
        coefTable[which(coefTable$Response %in% global$Response & coefTable$Predictor %in% global$Predictor), ]
      )
    
      i[, ncol(i)] <- ifelse(i$Response %in% global$Response & i$Predictor %in% global$Predictor, "c", "")
      
      i[, ncol(i) + 1] <- isSig(i$P.Value)
      
      names(i)[(ncol(i) - 1):ncol(i)] <- ""
      
      return(i) 
      
    } )
    
  }

  ret <- list(
    name = name,
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

  cat("\nStructural Equation Model of", x$name, "\n")
  
  cat("\nGroups =", x$group, ":", paste(names(x$group.coefs), collapse = ", "))
  
  cat("\n")
  
  cat("\nModel-wide Interactions:\n")
  
  cat("\n", captureTable(x$anovaInts))
  
  if(nrow(x$global) == 0) cat("\n No paths constrained to the global model") else {
    
    for(i in 1:nrow(x$global)) cat("\n", paste(x$global[i, "Predictor"], "->", x$global[i, "Response"], "constrained to the global model"), "\n")
  
  }
    
  cat("\n---\n\n")
  
  if(length(x$group.coefs) == 1) {
    
    cat("Global coefficients:\n")
    
    cat("\n", captureTable(x$group.coefs[[1]]), "\n")
    
    } else { 
      
      for(i in names(x$group.coefs)) { 
        
        cat(paste("Group", i), "coefficients: \n")
    
        cat("\n", captureTable(x$group.coefs[[i]]), "\n")
        
      }
  }
  
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05    c = constrained")
  
  cat("\n\n---\n\n")

  invisible(x)
  
}

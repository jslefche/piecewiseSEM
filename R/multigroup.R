#' Multigroup Analysis for Piecewise SEM
#' 
#' @param modelList a list of structural equations
#' @param group the name of the grouping variable in quotes
#' @param standardize The type of standardization: \code{none}, \code{scale}, \code{range}.
#' Default is \code{scale}.
#' @param standardize.type The type of standardized for non-Gaussian responses:
#' \code{latent.linear}, \code{Menard.OE}. Default is \code{latent.linear}.
#' @param test.type what kind of ANOVA should be reported. Default is type III
#' 
#' @author Jon Lefcheck <lefcheckj@@si.edu>
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
multigroup <- function(modelList, group, standardize = "scale", standardize.type = "latent.linear", test.type = "III") {
  
  name <- deparse(match.call()$modelList)
  
  data <- modelList$data
  
  modelList <- removeData(modelList, formulas = 1)
  
  intModelList <- lapply(modelList, function(i) {
    
    # rhs1 <- paste(paste(all.vars_trans(i)[-1], collapse = " + "))
    
    rhs2 <- paste(paste(all.vars_trans(i)[-1], "*", group), collapse = " + ")

    i <- update(i, formula(paste(". ~ ", rhs2))) # paste(rhs1, " + ", rhs2))))
    
    return(i)
    
  } )
  
  newModelList <- lapply(unique(data[, group]), function(i) update(as.psem(modelList), data = data[data[, group] == i, ]) )
  
  names(newModelList) <- unique(data[, group])
  
  coefsList <- lapply(newModelList, coefs, standardize, standardize.type, test.type)
  
  names(coefsList) <- unique(data[, group])
  
  coefTable <- coefs(modelList, standardize, standardize.type, test.type)

  anovaTable <- anova(as.psem(intModelList))[[1]]
  
  anovaInts <- anovaTable[grepl(":", anovaTable$Predictor), ]
  
  global <- anovaInts[anovaInts$P.Value >= 0.05, c("Response", "Predictor")]
  
  global$Predictor <- sub(":", "\\1", sub(group, "\\1", global$Predictor))
  
  if(nrow(global) == nrow(anovaInts)) newCoefsList <- list(global = coefTable) else {
  
    newCoefsList <- lapply(names(coefsList), function(i) {
      
      ct <- as.matrix(coefsList[[i]])

      idx <- which(apply(ct[, 1:2], 1, paste, collapse = "___") %in% apply(global[, 1:2], 1, paste, collapse = "___"))
      
      ct[idx, ] <- as.matrix(coefTable[idx, ])
      
      ct <- cbind(ct, ifelse(1:nrow(ct) %in% idx, "c", ""))
      
      # ct <- as.data.frame(ct)
      
      for(j in 1:nrow(ct)) {
        
        if(ct[j, ncol(ct)] == "c") {
      
          model <- modelList[[which(sapply(listFormula(modelList), function(x) all.vars.merMod(x)[1] == ct[j, "Response"]))]]
          
          data. <- data[data[, group] == i, ]
          
          sd.x <- GetSDx(model, modelList, data., standardize) 
          
          sd.x <- sd.x[which(names(sd.x) == ct[j, "Predictor"])]
          
          sd.y <- GetSDy(model, data., standardize, standardize.type)
          
          new.coef <- as.numeric(ct[j, "Estimate"]) * (sd.x/sd.y)
          
          ct[j, "Std.Estimate"] <- ifelse(length(new.coef) > 0, round(as.numeric(new.coef), 4), "-")
          
        }
        
      }
      
      ct <- as.data.frame(ct)
      
      ct[is.na(ct)] <- "-"  
    
      names(ct)[(ncol(ct) - 1):ncol(ct)] <- ""
      
      return(ct) 
      
    } )
    
    names(newCoefsList) <- names(coefsList)
    
  }
  
  if(nrow(global) == nrow(anovaInts)) gof <- fisherC(modelList) else {
    
    b <- basisSet(modelList)
    
    cf <- coefTable[coefTable$Response %in% global$Response & coefTable$Predictor %in% global$Predictor, ]
    
    b <- lapply(b, function(i) {
      
      for(j in 3:length(i)) {
        
        value <- cf[cf$Response == i[2] & cf$Predictor == i[j], "Estimate"]
        
        if(length(value) != 0) i[j] <- paste0("offset(", value, "*", i[j], ")")
        
      }
      
      return(i)
      
    } )
    
    if(length(b) == 0) b <- NULL
        
    gof <- fisherC(modelList, basis.set = b)
    
  }

  ret <- list(
    name = name,
    group = group,
    global = global,
    anovaInts = anovaInts,
    group.coefs = newCoefsList,
    Cstat = gof
  )
  
  class(ret) <- "multigroup.psem"
  
  return(ret)
  
}

#' Print multigroup
#' 
#' @param x an object to print
#' @param ... additional arguments to print
#' 
#' @method print multigroup.psem
#' 
#' @export
#' 
print.multigroup.psem <- function(x, ...) {

  cat("\nStructural Equation Model of", x$name, "\n")
  
  cat("\nGroups =", x$group, "[", paste(names(x$group.coefs), collapse = ", "), "]")
  
  cat("\n\n---\n")
  
  cat("\nGlobal goodness-of-fit:\n\n  Fisher's C =", as.character(x$Cstat[1]),
      "with P-value =", as.character(x$Cstat[3]),
      "and on", as.character(x$Cstat[2]), "degrees of freedom")
  
  cat("\n\n---\n")
  
  cat("\nModel-wide Interactions:\n")
  
  cat("\n", captureTable(x$anovaInts))
  
  if(nrow(x$global) == 0) cat("\n No paths constrained to the global model (P > 0.05)") else {
    
    for(i in 1:nrow(x$global)) cat("\n", paste(x$global[i, "Predictor"], "->", x$global[i, "Response"], "constrained to the global model"))
  
  }
  
  cat("\n\n---\n\n")
  
  if(length(x$group.coefs) == 1) {
    
    cat("Global coefficients:\n")
    
    cat("\n", captureTable(x$group.coefs[[1]]), "\n")
    
    } else { 
      
      for(i in names(x$group.coefs)) { 
        
        cat(paste0("Group [", i, "]"), "coefficients: \n")
    
        cat("\n", captureTable(x$group.coefs[[i]]), "\n")
        
      }
  }
  
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05    c = constrained")

  # if(length(x$LRTs) > 0) {
  #   
  #   cat("\n\n---\n\n")
  #   
  #   cat("Chi-squared difference tests:\n")
  #   
  #   for(i in names(x$LRTs)) {
  #     
  #     cat(paste("\n  Global model vs. Group", i))
  #     
  #     cat("\n  ", captureTable(x$LRTs[[i]][[i]]), "\n")
  #         
  #   }
  #   
  # }

  invisible(x)
  
}

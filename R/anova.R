#' Compare models using Fisher's C scores

#' @param object a \code{\link{psem}} object
#' @param object2 an optional second \code{\link{psem}} object
#' @param fun what anova function to use?
#' @param ... options for anova functions 

anova.psem <- function(object, object2 = NULL, fun = Anova, ...) {
  require(car)
#  dots <- list(...)
  
   #are we getting the LRT/F tables from
   if(is.null(object2) ){
      object <- removeData(object, formulas = 1)
      
      #get residuals of relationships
      ret <- lapply(object, fun, options)
      
      #get column names
      formulaList <- listFormula(object)
      names(ret) <- sapply(formulaList, function(x) all.vars(x)[1])
      
      ret

   }else{
      anova.psemlist(object, object2)
    }   
}


anova.psemlist <- function(mod1, mod2){
  ret <- rbind(fisherC(mod1), fisherC(mod2))
  ret <- rbind(ret, abs(ret[1,] - ret[2,]))
  
  ret[3,3] <- 1-pchisq(ret[3,1], df = ret[3,2])
  
  rownames(ret) <- c(1, 2, "Difference")
  
  ret
  
}
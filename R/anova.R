#' Compare models using Fisher's C scores

#' @param object a \code{\link{psem}} object
#' @param object2 an optional second \code{\link{psem}} object
#' @param fun what anova function to use. Default is car::Anova
#' @param ... options for anova functions
#' 
#' @export
#'
anova.psem <- function(object, object2 = NULL, fun = Anova, ...) {

  # require(car)
  # dots <- list(...)

  #are we getting the LRT/F tables from
  if(is.null(object2) ){

   object <- removeData(object, formulas = 1)

   #get residuals of relationships
   ret <- lapply(object, fun, ...)

    #get column names
    formulaList <- listFormula(object)

    names(ret) <- sapply(formulaList, function(x) all.vars(x)[1])

    } else {

      ret <- rbind(fisherC(object), fisherC(object2))
      
      ret <- rbind(ret, abs(ret[1, ] - ret[2, ]))
      
      ret[3, 3] <- 1 - pchisq(ret[3,1], df = ret[3, 2])
      
      rownames(ret) <- c(1, 2, "Difference")

    }
  
  return(ret)

}

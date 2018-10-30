#' Chi-squared Difference Test for Model Comparison
#' 
#' @param object a \code{\link{psem}} object
#' @param object2 an optional second \code{\link{psem}} object
#' @param fun what anova function to use. Default is car::Anova
#' @param ... additional arguments for car::Anova
#' 
#' @return a table reporting the Fisher's C statistics for each model, and the result of 
#' the Chi-squared difference test
#' 
#' @method anova psem
#' 
#' @export
#'
anova.psem <- function(object, object2 = NULL, ...) {

  if(is.null(object2)) {

   object <- removeData(object, formulas = 1)

   #get residuals of relationships
   ret <- lapply(object, car::Anova, ...)
   
   ret <- do.call(rbind, lapply(ret, function(i) {
     
     response <- gsub("Response: ", "", attr(i, "heading")[grepl("Response:", attr(i, "heading"))])
     
     dat <- as.data.frame(i)
     
     predictors <- rownames(dat)[-nrow(dat)]
     
     data.frame(
       Response = response,
       Predictor = predictors,
       Test.stat = dat[-nrow(dat), 1],
       DF = dat[-nrow(dat), "Df"],
       P.Value = dat[-nrow(dat), ncol(dat)]
     )
     
   } ) )
   
  } else {
    
    ret <- rbind(fisherC(object), fisherC(object2))
    
    ret <- rbind(ret, abs(ret[1, ] - ret[2, ]))
    
    ret[3, 3] <- 1 - pchisq(ret[3, 1], df = ret[3, 2])
    
    rownames(ret) <- c(1, 2, "Difference")
    
    }
  
  return(ret)

}

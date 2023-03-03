#' ANOVA and chi-squared difference test for model comparison
#' 
#' Compute analysis of variance table for one or more structural equation models.
#' 
#' Additional models will be tested against the first model using a Chi-squared difference test.
#' 
#' @param object a \code{psem} object
#' @param ... additional objects of the same type
#' @param anovafun The function used for ANOVA. Defaults to \code{Anova}
#' @param digits number of digits to round results. Default is 3
#' 
#' @return an F, LRT, or other table for a single model, or a list of
#'  comparisons between multiple models
#' 
#' @author Jon Lefcheck <lefcheckj@@si.edu>, Jarrett Byrnes <jarrett.byrnes@@umb.edu>  
#' 
#' @seealso \code{\link{Anova}}
#' 
#' @examples 
#' data(keeley)
#' 
#' mod1 <- psem(
#' lm(rich ~ cover, data = keeley),
#' lm(cover ~ firesev, data = keeley),
#' lm(firesev ~ age, data = keeley),
#' data = keeley
#' )
#' 
#' # get type II Anova
#' anova(mod1)
#' 
#' # conduct LRT
#' mod2 <- psem(
#'   lm(rich ~ cover, data = keeley),
#'   lm(cover ~ firesev, data = keeley),
#'   age ~ 1,
#'   data = keeley
#' )
#' 
#' anova(mod1, mod2)
#' 
#' @method anova psem
#' 
#' @export
#' 
anova.psem <- function(object, ..., digits = 3, anovafun = "Anova") {
 
  nm <- deparse(substitute(object))
  
  nms <- deparse(substitute(...))
  
  dots <- list(object, ...)
  
  if(length(dots) > 1) names(dots) <- c(nm, nms) else names(dots) <- nm
  
  if(length(dots) > 1) {
     
    anovaLRT(dots)
   
    } else {
     
      anovaTable(object, anovafun = anovafun, digits = digits)
       
       }
}

#' Single anova
#' 
#' @keywords internal
#' 
anovaTable <- function(object, anovafun = "Anova", digits = 3) {
  
  object <- removeData(object, formulas = 1)
  
  if(anovafun == "Anova") a <- car::Anova else stop("Unsupported ANOVA function")
    
    # if(anovafun == "aov") a <- aov else 
  
  tests <- lapply(object, function(x) a(x))
  
  names(tests) <- get_response(object)
  
  ret <- list(do.call(rbind, lapply(names(tests), function(x) {
    
    i <- tests[[x]]
    
    response <- x
    
    dat <- as.data.frame(i)
    
    predictors <- rownames(dat)
    
    DF <- ifelse(any(grepl("numDF", colnames(dat))), dat$numDF, dat$Df)
    
    Test.Stat <- ifelse(any(grepl("F-value", colnames(dat))), dat$`F-value`, dat[, 1])
    
    ret <- data.frame(
      Response = response,
      Predictor = predictors,
      Test.Stat = round(Test.Stat, 1),
      DF = DF,
      P.Value = round(dat[, ncol(dat)], 4)
    )
    
    ret <- ret[!ret$Predictor %in% c("(Intercept)", "Residuals"), ]
    
    ret <- cbind.data.frame(ret, isSig(ret$P.Value))
    
    colnames(ret)[ncol(ret)] <- ""
    
    rownames(ret) <- NULL
    # ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)
    
    return(ret)
    
  } ) ) )
  
  class(ret) <- "anova.psem"
  
  return(ret)
  
}

#' Chi-square difference test
#' 
#' @keywords internal
#' 
anovaLRT <- function(object) {

  model1 <- object[[1]]
  
  ChiSq1 <- LLchisq(model1)
  
  ret1 <- data.frame(
    Df = ChiSq1$df,
    AIC(model1),
    Chisq = ChiSq1$Chisq,
    Chisq.diff = NA,
    Df.diff = NA,
    P.value = NA,
    sig = NA
  )
  
  # rownames(ret1) <- deparse(substitute(object[[1]]))
  
  ret2 <- do.call(rbind, lapply(2:length(object), function(i) {
    
    model2 <- object[[i]]
  
    ChiSq2 <- LLchisq(model2)
    
    Chisq.diff <- abs(ChiSq1$Chisq - ChiSq2$Chisq)
    
    df.diff <- abs(ChiSq1$df - ChiSq2$df)
    
    pvalue <- round(1 - pchisq(Chisq.diff, df = df.diff), 4)
    
    ret2 <- data.frame(
      Df = ChiSq2$df,
      AIC(model2),
      Chisq = ChiSq2$Chisq,
      Chisq.diff = Chisq.diff,
      Df.diff = df.diff,
      P.value = pvalue,
      sig = isSig(pvalue)
    )
    
    # rownames(ret2) <- deparse(substitute(object[[i]]))
    
    return(ret2)
    
  } ) )

  ret <- rbind(ret1, ret2)
  
  colnames(ret)[ncol(ret)] <- ""
  
  rownames(ret) <- c(names(object)[1], paste("vs", names(object)[-1]))
  
  ret[is.na(ret)] <- ""
  
  ret <- list(ret)
  
  class(ret) <- "anova.psem"
  
  return(ret)

}

#' Print anova
#' 
#' @param x an object of class anova.psem
#' @param ... further arguments passed to or from other methods
#' 
#' @method print anova.psem
#' 
#' @export
#' 
print.anova.psem <- function(x, ...) {
  
  if(any(grepl("Response", names(x[[1]])))) {
    
    cat(captureTable(x[[1]]))
    
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05") 
    
    } else {
    
    cat("Chi-square Difference Test\n")
      
    cat("\n")
    
    cat(captureTable(x[[1]], row.names = TRUE))
    
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05")
    
    }
  
}

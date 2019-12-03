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
#' @author Jon Lefcheck <LefcheckJ@@si.edu>, Jarrett Byrnes <jarrett.byrnes@@umb.edu>  
#' 
#' @seealso \code{\link{car::Anova}}
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
 
  dots <- list(object, ...)
  
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
  
  model1.summary <- fisherC(model1, .progressBar = FALSE)
  
  ret1 <- data.frame(
    AIC = AIC(model1),
    BIC = BIC(model1),
    Fisher.C = model1.summary$Fisher.C,
    Fisher.C.Diff = NA,
    DF.diff = NA,
    P.value = NA,
    sig = NA
  )
  
  colnames(ret1)[ncol(ret1)] <- ""
  
  rownames(ret1) <- "1"
  
  ret2 <- do.call(rbind, lapply(2:length(object), function(i) {
    
    model2 <- object[[i]]
  
    model2.summary <- fisherC(model2, .progressBar = FALSE)
    
    Cdiff <- abs(model1.summary$Fisher.C - model2.summary$Fisher.C)
    
    dfdiff <- abs(model1.summary$df - model2.summary$df)
    
    pvalue <- round(1 - pchisq(Cdiff, df = dfdiff), 4)
    
    ret <- data.frame(
      AIC = AIC(model2),
      BIC = BIC(model2),
      Fisher.C = model2.summary$Fisher.C,
      Fisher.C.Diff = Cdiff,
      DF.diff = dfdiff,
      P.value = pvalue,
      sig = isSig(pvalue)
    )
     
    colnames(ret)[ncol(ret)] <- ""
    
    rownames(ret) <- c(paste("vs", i))
    
    return(ret)
    
  } ) )

  ret <- rbind(ret1, ret2)
  
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

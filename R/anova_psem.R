#' ANOVA and chi-squared difference test for model comparison
#' 
#' Compute analysis of variance table for one or more structural equation models.
#' 
#' @param mod a \code{psem} object
#' @param mod2 a \code{psem} object for comparison. Defaults to NULL to allow for a LRT or F tables of all fit model pieces.
#' @param anovafun The function used for ANOVA. Defaults to \code{\link{car::Anova}}
#' @param digits number of digits to round results. Default is 3
#' @param ... additional arguments passed to \code{anovafun}
#' 
#' @return an F, LRT, or other table for a single model, or a list of
#'  comparisons between multiple models
#' 
#' @author Jon Lefcheck <LefcheckJ@@si.edu>, Jarrett Byrnes <jarrett.byrnes@@umb.edu>  
#' 
#' @seealso The model fitting function \code{\link{psem}}
#' 
#' @examples 
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
anova.psem <- function(mod, mod2 = NULL, digits = 3, anovafun = car::Anova, ...) {
 
   if(!is.null(mod2)) {
     
    anovaLRT(mod, mod2)
     
     } else {
       
      anovaTable(mod, anovafun = anovafun, digits = digits, ...)
       
       }
}

#' Single anova
#' 
#' @keywords internal
#' 
anovaTable <- function(object, anovafun = car::Anova, digits = 3, ...) {
  
  object <- removeData(object, formulas = 1)
  
  tests <- lapply(object, function(x) anovafun(x, ...))
  
  names(tests) <- get_response(object)
  
  ret <- list(do.call(rbind, lapply(names(tests), function(x) {
    
    i <- tests[[x]]
    
    response <- x
    
    dat <- as.data.frame(i)
    
    predictors <- rownames(dat)[-nrow(dat)]
    
    ret <- data.frame(
      Response = response,
      Predictor = predictors,
      Test.Stat = round(dat[-nrow(dat), 1], 1),
      DF = dat[-nrow(dat), "Df"],
      P.Value = round(dat[-nrow(dat), ncol(dat)], 4)
    )
    
    ret <- cbind.data.frame(ret, isSig(ret$P.Value))
    
    colnames(ret)[ncol(ret)] <- ""
    
    rownames(ret) <- NULL
    # ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)
    
    return(ret)
    
  } ) ) )
  
  class(ret) <- "anova.psem"
  
  return(ret)
  
}

#' Chi-square differnce test
#' 
#' @keywords internal
#' 
anovaLRT <- function(...) {

  dots <- list(...)
  
  combos <- combn(length(dots), 2)
    
  combos <- split(t(combos), seq(ncol(combos)))
  
  ret <- lapply(combos, function(i) {
    
    nm1 <- deparse(substitute(dots[[i[1]]]))
    
    nm2 <- deparse(substitute(dots[[i[2]]]))
    
    model1 <- dots[[i[1]]]
    
    model2 <- dots[[i[2]]]
  
    model1.summary <- fisherC(model1, .progressBar = FALSE)
    
    model2.summary <- fisherC(model2, .progressBar = FALSE)
    
    Cdiff <- abs(model1.summary$Fisher.C - model2.summary$Fisher.C)
    
    dfdiff <- abs(model1.summary$df - model2.summary$df)
    
    pvalue <- round(1 - pchisq(Cdiff, df = dfdiff), 4)
    
    ret <- data.frame(
      AIC = c(AIC(model1), AIC(model2)),
      BIC = c(BIC(model1), BIC(model2)),
      # BIC = c(model1.summary$IC$BIC, model2.summary$IC$BIC),
      Fisher.C = c(model1.summary$Fisher.C, model2.summary$Fisher.C),
      Fisher.C.Diff = c("", Cdiff),
      DF.diff = c("", dfdiff),
      P.value = c("", pvalue),
      sig = c("", isSig(pvalue))
    )
     
    colnames(ret)[ncol(ret)] <- ""
    
    rownames(ret) <- c(paste(i[1]), paste("vs", i[2]))
    
    return(ret)
    
  } )

  class(ret) <- "anova.psem"
  
  return(ret)

}

#' Print anova
#' 
#' @param x an object of class anova.psem
#' 
#' @method print anova.psem
#' 
#' @export
#' 
print.anova.psem <- function(x) {
  
  if(grepl("Response", colnames(x[[1]])[1])) print(x[[1]]) else {
  
    cat("Chi-square Difference Test\n")
    
    cat("\n")
    
    lapply(x, function(i) { print(i); cat("\n") })
    
  }
  
}

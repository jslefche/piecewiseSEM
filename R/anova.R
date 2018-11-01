#' ANOVA and chi-squared difference test for model comparison
#' 
#' Compute analysis of variance table for one or more structural equation models.
#' 
#' @param object a \code{psem} object
#' @param ... additional \code{psem} objects
#' @param test.type what kind of ANOVA should be reported. Default is type II
#' 
#' @return an ANOVA table for a single model, a list of comparisons between multiple
#' models
#' @author Jon Lefcheck <LefcheckJ@@si.edu>, Jarrett Byrnes
#' 
#' @seealso The model fitting function \code{\link{psem}}
#' 
#' @method anova psem
#' 
#' @export
#' 
anova.psem <- function(object, ..., 
                       test.type = "II") {

  l <- list(object, ...)
  
  if(length(l) == 1) {

   object <- removeData(object, formulas = 1)

   ret <- lapply(object, function(x) car::Anova(x, type = test.type))
   
   ret <- list(do.call(rbind, lapply(ret, function(i) {
     
     response <- gsub("Response: ", "", attr(i, "heading")[grepl("Response:", attr(i, "heading"))])
     
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
     
     names(ret)[ncol(ret)] <- ""
     
     # ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)
     
     return(ret)
     
   } ) ) )
   
  } else {
    
    # need to stop if fit to different datasets
    # lapply(l, function(x) GetData(x))
    
    combos <- combn(length(l), 2)
    
    combos <- split(t(combos), seq(ncol(combos)))
    
    ret <- lapply(combos, function(i) {
      
      nm1 <- deparse(substitute(l[[i[1]]]))
      
      nm2 <- deparse(substitute(l[[i[2]]]))
      
      model1 <- l[[i[1]]]
      
      model2 <- l[[i[2]]]
    
      model1.summary <- summary(model1, .progressBar = FALSE)
      
      model2.summary <- summary(model2, .progressBar = FALSE)
      
      Cdiff <- abs(model1.summary$Cstat$Fisher.C - model2.summary$Cstat$Fisher.C)
      
      dfdiff <- abs(model1.summary$Cstat$df - model2.summary$Cstat$df)
      
      pvalue <- 1 - pchisq(Cdiff, df = dfdiff)
      
      ret <- data.frame(
        AIC = c(model1.summary$IC$AIC, model2.summary$IC$AIC),
        BIC = c(model1.summary$IC$BIC, model2.summary$IC$BIC),
        Fisher.C = c(model1.summary$Cstat$Fisher.C, model2.summary$Cstat$Fisher.C),
        Fisher.C.Diff = c("", Cdiff),
        DF.diff = c("", dfdiff),
        P.value = c("", pvalue),
        sig = c("", isSig(pvalue))
      )
       
      colnames(ret)[ncol(ret)] <- ""
      
      rownames(ret) <- c(paste(i[1]), paste("vs", i[2]))
      
      return(ret)
      
    } )
    
  }

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
  
  if(length(x) == 1) print(x[[1]]) else {
  
    cat("Chi square difference test\n")
    
    cat("\n")
    
    lapply(x, function(i) { print(i); cat("\n") })
    
  }
  
}

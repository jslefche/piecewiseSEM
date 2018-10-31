#' Chi-squared Difference Test for Model Comparison
#' 
#' @param object a \code{psem} object
#' @param ... additional \code{psem} objects
#' @param test.type 
#' 
#' @return 
#' 
#' @method anova psem
#' 
#' @export
#'
anova.psem <- function(object, ..., test.type = "II") {

  l <- list(object, ...)
  
  if(length(l) == 1) {

   object <- removeData(object, formulas = 1)

   ret <- lapply(object, function(x) car::Anova(x, type = test.type))
   
   ret <- do.call(rbind, lapply(ret, function(i) {
     
     response <- gsub("Response: ", "", attr(i, "heading")[grepl("Response:", attr(i, "heading"))])
     
     dat <- as.data.frame(i)
     
     predictors <- rownames(dat)[-nrow(dat)]
     
     data.frame(
       Response = response,
       Predictor = predictors,
       Test.Stat = dat[-nrow(dat), 1],
       DF = dat[-nrow(dat), "Df"],
       P.Value = dat[-nrow(dat), ncol(dat)]
     )
     
   } ) )
   
  } else {

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
  
  return(ret)

}

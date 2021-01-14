#' Generalized chi-squared for piecewise SEM
#' 
#' Derivation of log-likelihoods to be used in determining the
#' goodness-of-fit for piecewise structural equation models.
#' 
#' Here, a list of saturated models is first derived from the list of
#' structured equations using the basis set. Then, the differences in summed
#' log-likelihoods is computed and used to calculate  the Chi-squared statistic.
#' 
#' @param modelList A list of structural equations created using \code{psem}.
#' @param basis.set An optional list of independence claims.
#' @param interactions whether interactions should be included in basis set. 
#' Default is FALSE
#' 
#' @return a data.frame corresponding to the Chi-squared statistic, d.f., and P-value
#' 
#' @author Jon Lefcheck <LefcheckJ@@si.edu>
#' 
#' @seealso \code{\link{basisSet}}
#' 
#' @references Shipley, Bill, and Jacob C. Douma. "Generalized AIC and chi‐squared statistics 
#' for path models consistent with directed acyclic graphs." Ecology 101.3 (2020): e02960.
#' 
#' @examples 
#' mod <- psem(
#' lm(rich ~ cover, data = keeley),
#' lm(cover ~ firesev, data = keeley),
#' lm(firesev ~ age, data = keeley),
#' data = keeley
#' )
#' 
#' LLchisq(mod)
#' 
#' @export 
#' 
LLchisq <- function(modelList, basis.set = NULL, interactions = FALSE) {
  
  if(is.null(basis.set)) b <- basisSet(modelList, direction = NULL, interactions) else b <- basis.set

  if(length(b) == 0) {
    
    data.frame(Chisq = 0, df = 0, P.Value = 1)
    
  } else {
  
    if(class(modelList) == "psem") data <- modelList$data else data <- GetData(modelList)
    
    modelList <- removeData(modelList, formulas = 1)
  
    satModelList <- getSatModels(b, modelList, data)
    
    M1 <- sapply(modelList, logLik)
    
    M2 <- sapply(satModelList, logLik)
    
    ChiSq_ML <- -2*sum(M1 - M2)
    
    if(ChiSq_ML < 0) {
      
      ChiSq_ML <- NA
      
      warning("Check convergence: log-likelihood estimates lead to negative Chi-squared!")
      
    }
    
    DF <- sum(sapply(satModelList, function(x) attributes(logLik(x))$df)) -
      sum(sapply(modelList, function(x) attributes(logLik(x))$df))
    
    P <- 1 - pchisq(ChiSq_ML, DF)
    
    ret <- data.frame(Chisq = ChiSq_ML, df = DF, P.Value = P)
    
    ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 3)
    
    return(ret)
  
  }
  
}

getSatModels <- function(b, modelList, data) {
  
  lapply(1:length(modelList), function(i) {
    
    model <- modelList[[i]]
    
    y <- all.vars_trans(model)[1]
    
    newVars <- na.omit(unlist(sapply(b, function(x) if(x[2] == y) x[-c(2)] else NA)))
    
    newVars <- newVars[!duplicated(newVars)]
    
    newVars <- newVars[!newVars %in% all.vars_trans(model)]
    
    if(length(newVars) == 0) satModel <- model else {
      
      satModel <- suppressWarnings(
        update(model,
               formula(paste(". ~ . +", paste(newVars, collapse = " + "))),
               data = data) )
      
    }
    
  } )
  
}
      
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
#' @export 
#' 
LLchisq <- function(modelList, basis.set = NULL, interactions = FALSE) {
  
  if(is.null(basis.set)) b <- basisSet(modelList, direction = NULL, interactions) else b <- basis.set
  
  if(any(duplicated(names(b))) & conserve == FALSE & is.null(direction)) dupOutput(b)
  
  if(length(b) == 0) {
    
    data.frame()
    
  } else {
  
    if(class(modelList) == "psem") data <- modelList$data else data <- GetData(modelList)
    
    modelList <- removeData(modelList, formulas = 1)
  
    # loop over modelList and update models to get saturated modelw
    satModelList <- getSatModels(b, modelList, data)
    
    # get LL from model List
    
    M1 <- sapply(modelList, logLik)
    
    M2 <- sapply(satModelList, logLik)
    
    ChiSq <- sum(-2*(M1 - M2))
    
    DF <- sum(sapply(satModelList, function(x) attributes(logLik(x))$df)) -
      sum(sapply(modelList, function(x) attributes(logLik(x))$df))
    
    P <- 1 - pchisq(ChiSq, DF)
    
    ret <- data.frame(Chisq = ChiSq, df = DF, P.Value = P)
    
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
      
      if(any(class(model) %in% c("lmerMod", "merModLmerTest", "lmerModLmerTest", "glmerMod"))) {
        
        satModel <- suppressWarnings(
          update(model,
                 formula(paste(". ~ . +", paste(newVars, collapse = " + "), " + ", onlyBars(formula(bMod)))),
                 data = data)
        )
        
      } else {
        
        satModel <- suppressWarnings(
          update(model,
                 formula(paste(". ~ . +", paste(newVars, collapse = " + "))),
                 data = data)
        )
        
      }
      
    }
    
  } )
  
}
      
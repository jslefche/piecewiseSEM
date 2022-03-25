#' Information criterion values for SEM
#' 
#' @param modelList a list of structural equations
#' @param AIC.type whether the log-likelihood \code{"loglik"} or d-sep \code{"dsep"} AIC score 
#' should be reported. Default is \code{"loglik"}
#' @param Cstat Fisher's C statistic obtained from \code{fisherC}
#' @param add.claims an optional vector of additional independence claims (P-values) 
#' to be added to the basis set
#' @param basis.set An optional list of independence claims.
#' @param direction a vector of claims defining the specific directionality of any independence 
#' claim(s)
#' @param interactions whether interactions should be included in independence claims. 
#' Default is FALSE
#' @param conserve whether the most conservative P-value should be returned (See Details) 
#' Default is FALSE
#' @param conditional whether the conditioning variables should be shown in the table. 
#' Default is FALSE
#' @param .progressBar an optional progress bar. Default is FALSE
#' 
#' @return a data.frame of AIC, AICc, d.f., and sample size
#' 
#' @author Jon Lefcheck <LefcheckJ@@si.edu>, Jim Grace
#' 
#' @references Shipley, Bill, and Jacob C. Douma. "Generalized AIC and chi‐squared statistics 
#' for path models consistent with directed acyclic graphs." Ecology 101.3 (2020): e02960.
#' 
#' Shipley, Bill. "The AIC model selection method applied to path analytic models compared using 
#' a d‐separation test." Ecology 94.3 (2013): 560-564.
#' 
#' @export
#' 
AIC_psem <- function(modelList, AIC.type = "loglik", 
                     Cstat = NULL, add.claims = NULL, basis.set = NULL, direction = NULL, 
                     interactions = FALSE, conserve = FALSE, conditional = FALSE, 
                     .progressBar = FALSE) {
  
  modelList <- removeData(modelList, formulas = 1)
  
  K <- do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))
  
  n.obs <- min(sapply(modelList, nObs))
  
  if(AIC.type == "loglik") {
    
    aic <- sum(sapply(modelList, AIC))
    
    aicc <- sum(sapply(modelList, MuMIn::AICc))
    
  } else if(AIC.type == "dsep") {
    
    if(missing(Cstat)) Cstat <- fisherC(modelList, add.claims, basis.set, direction, conserve, conditional, .progressBar)
    
    aic <- as.numeric(Cstat[1] + 2*K)
    
    aicc <- as.numeric(Cstat[1]) + 2*K*(n.obs/(n.obs - K - 1))
    
    # BIC <- as.numeric(Cstat[1] + log(n.obs)*K)
    
  }
  
  ret <- data.frame(
    AIC = aic,
    AICc = aicc,
    # BIC = bic,
    K = K,
    n = n.obs
  )
  
  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 3)
  
  return(ret)
  
}

#' Generic function for SEM AIC(c) score
#'
#' @param object a psem object
#' @param ... additional arguments to AIC
#' @param AIC.type whether the log-likelihood \code{"loglik"} or d-sep \code{"dsep"} AIC score 
#' should be reported. Default is \code{"loglik"}
#' @param aicc whether correction for small sample size should be applied. Default is \code{FALSE}
#'  
#' @method AIC psem
#' 
#' @examples 
#' 
#' mod <- psem(
#' lm(rich ~ cover, data = keeley),
#' lm(cover ~ firesev, data = keeley),
#' lm(firesev ~ age, data = keeley),
#' data = keeley
#' )
#' 
#' # Get log-likelihood based AIC
#' AIC(mod)
#' 
#' # Get d-sep based AIC
#' AIC(mod, AIC.type = "dsep")
#'  
#' @export
#' 
AIC.psem <- function(object, ..., AIC.type = "loglik", aicc = FALSE) {
  
  aicx <- AIC_psem(object, AIC.type)
  
  if(missing(...)) ret <- aicx else {
    
    dots <- list(object, ...)
    
    aicy <- do.call(rbind, lapply(2:length(dots), function(i) AICy <- AIC_psem(dots[[i]], AIC.type)))

    ret <- rbind(aicx, aicy)
    
  }
  
  if(aicc == FALSE) ret <- subset(ret, select = -c(AICc))

  return(ret)
  
}

#' Information criterion values for SEM
#' 
#' @param modelList a list of structural equations
#' @param Cstat Fisher's C statistic obtained from \code{fisherC}
#' @param add.claims an optional vector of additional independence claims (P-values) 
#' to be added to the basis set
#' @param direction a vector of claims defining the specific directionality of any independence 
#' claim(s)
#' @param conserve whether the most conservative P-value should be returned (See Details) 
#' Default is FALSE
#' @param conditional whether the conditioning variables should be shown in the table. 
#' Default is FALSE
#' @param .progressBar an optional progress bar. Default is FALSE
#' 
#' @return a vector of AIC, AICc, BIC, d.f., and sample size
#' 
infCrit <- function(modelList, Cstat, add.claims = NULL, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = FALSE) {

  if(missing(Cstat)) Cstat <- fisherC(modelList, add.claims, direction, conserve, conditional, .progressBar)

  modelList <- removeData(modelList, formulas = 1)

  K <- do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))

  n.obs <- min(sapply(modelList, nObs))

  pwAIC <- as.numeric(Cstat[1] + 2*K)

  pwAICc <- pwAIC * (n.obs/(n.obs - K - 1))

  pwBIC <- as.numeric(Cstat[1] + log(n.obs)*K)

  ret <- data.frame(
    AIC = pwAIC,
    AICc = pwAICc,
    BIC = pwBIC,
    K = K,
    n = n.obs
  )

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 3)

  return(ret)

}

#' Generalized function for SEM AIC(c) score
#'
#' @param object a psem object
#' @param ... additional arguments to AIC
#' @param aicc whether correction for small sample size should be applied. Default is \code{FALSE}
#'  
#' @method AIC psem
#'  
#' @export
#' 
AIC.psem <- function(object, ..., aicc = FALSE) {

  aicx <- infCrit(object)

  if(aicc == FALSE) AICx <- aicx$AIC else AICx <- aicx$AICc

  if(missing(...)) ret <- AICx else {

    aicy <- infCrit(...)

    if(aicc == FALSE) AICy <- aicy$AIC else AICy <- aicy$AICc

    dfx <- aicx$K; dfy <- aicy$K

    ret <- data.frame(
      df = c(dfx, dfy),
      AIC = c(AICx, AICy),
      row.names = c(deparse(substitute(x)), deparse(substitute(y)))
    )

  }

  return(ret)

}

#' Generalized function for SEM BIC score
#'
#' @param object a psem object
#' @param ... additional arguments to BIC
#' 
#' @method BIC psem
#' 
#' @export
#' 
BIC.psem <- function(object, ...) infCrit(object)$BIC

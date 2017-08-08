#' Return information criterion values for piecewise SEM
#'
#' @param modelList a list of structural equations

infCrit <- function(modelList, Cstat = NULL, add.claims = NULL, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = FALSE) {

  if(is.null(Cstat)) Cstat <- fisherC(modelList, add.claims, direction, conserve, conditional, .progressBar)

  modelList <- removeData(modelList, formulas = 1)

  K <- do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))

  n.obs <- min(sapply(modelList, nobs.))

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

#` Generalized function for extraction AIC(c) score
AIC.psem <- function(modelList) infCrit(modelList)$AIC

AICc.psem <- function(modelList) infCrit(modelList)$AICc

#` Generalized function for extraction BIC score
BIC.psem <- function(modelList) infCrit(modelList)$BIC


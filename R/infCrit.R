#' Return information criterion values for piecewise SEM
#'
#' @param modelList a list of structural equations

infCrit <- function(modelList, Cstat, add.claims = NULL, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = FALSE) {

  if(missing(Cstat)) Cstat <- fisherC(modelList, add.claims, direction, conserve, conditional, .progressBar)

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

#' @export AIC.psem
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

#` Generalized function for extraction BIC score
BIC.psem <- function(object, ...) infCrit(object)$BIC

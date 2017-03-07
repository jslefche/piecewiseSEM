#' Summarize tests of directed separation using Fisher's C statistic

#' @param dTable a list of structural equations

InfCrit <- function(modelList, C) {

  mList <- modelList[!sapply(modelList, function(i) any(class(i) %in% c("formula.cerror")))]

  K <- do.call(sum, lapply(mList, function(i) attr(logLik(i), "df")))

  n.obs <- min(sapply(mList, nobs))

  pwAIC <- as.numeric(C[1] + 2*K)

  pwAICc <- pwAIC * (n.obs/(n.obs - K - 1))

  pwBIC <- as.numeric(C[1] + log(n.obs)*K)

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

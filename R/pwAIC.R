#' Summarize tests of directed separation using Fisher's C statistic

#' @param dTable a list of structural equations

pwAIC <- function(modelList, C) {

  mList <- modelList[!sapply(modelList, function(i) class(i) %in% c("formula.cerror"))]

  K <- do.call(sum, lapply(mList, function(i) attr(logLik(i), "df")))

  pwAIC <- as.numeric(C[1] + 2*K)

  n.obs <- min(sapply(mList, nobs))

  pwAICc <- pwAIC * (n.obs/(n.obs - K - 1))

  ret <- data.frame(
    AIC = pwAIC,
    AICc = pwAICc,
    K = K,
    n = n.obs
  )

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 3)

  return(ret)

}

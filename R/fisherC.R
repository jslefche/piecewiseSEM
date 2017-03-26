#' Summarize tests of directed separation using Fisher's C statistic
#'
#' @param dTable a list of structural equations

fisherC <- function(dTable) {

  if(length(dTable) == 0) {

    C <- 0

    p <- 1

    DF <- 0

  } else {

    if(any(dTable$P.value == 0)) dTable$P.value = dTable$P.Value + 1e-20

    C <- -2*sum(log(dTable$P.Value))

    p <- 1 - pchisq(C, 2*nrow(dTable))

    DF <- 2*nrow(dTable)

  }

  ret <- data.frame(Fisher.C = C, df = DF, P.Value = p)

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 3)

  return(ret)

}

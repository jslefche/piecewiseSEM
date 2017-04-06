#' Summarize tests of directed separation using Fisher's C statistic
#'
#' @param dTable a list of structural equations

fisherC <- function(dTable, claims = NULL) {

  if(length(dTable) == 0) {

    C <- 0

    DF <- 0

    P <- 1

  } else {

    ps <- dTable$P.Value

    if(!is.null(claims)) {

      ps <- ps + claims

      message("Fisher's C has been adjusted to include additional claims not shown in the tests of directed separation.")

      }

    if(any(ps == 0)) ps <- ps + 1e-20

    C <- -2 * sum(log(ps))

    DF <- 2 * length(ps)

    P <- 1 - pchisq(C, 2 * DF)

  }

  ret <- data.frame(Fisher.C = C, df = DF, P.Value = P)

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 3)

  return(ret)

}

#' Summarize tests of directed separation using Fisher's C statistic
#'
#' @param dTable a list of structural equations

fisherC <- function(dTable, add.claims = NULL, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = FALSE) {

  if(class(dTable) == "psem") dTable <- dSep(dTable, direction, conserve, conditional, .progressBar)

  if(length(dTable) == 0) {

    Cstat <- 0

    DF <- 0

    P <- 1

  } else {

    ps <- dTable$P.Value

    if(!is.null(add.claims)) {

      ps <- ps + add.claims

      message("Fisher's C has been adjusted to include additional claims not shown in the tests of directed separation.")

      }

    if(any(ps == 0)) ps <- ps + 1e-20

    Cstat <- -2 * sum(log(ps))

    DF <- 2 * length(ps)

    P <- 1 - pchisq(Cstat, DF)

  }

  ret <- data.frame(Fisher.C = Cstat, df = DF, P.Value = P)

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 3)

  return(ret)

}

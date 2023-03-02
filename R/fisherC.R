#' Summarize tests of directed separation using Fisher's C statistic
#'
#' @param dTable a \code{data.frame} containing tests of directed separation from \code{dSep}
#' @param add.claims an optional vector of additional independence claims (i.e., P-values) 
#' to be added to the basis set
#' @param basis.set An optional list of independence claims.
#' @param direction a vector of claims defining the specific directionality of any independence 
#' claim(s)
#' @param interactions whether interactions should be included in independence claims. 
#' Default is FALSE
#' @param conserve whether the most conservative P-value should be returned. 
#' Default is FALSE
#' @param conditional whether the conditioning variables should be shown in the table. 
#' Default is FALSE
#' @param .progressBar an optional progress bar. Default is FALSE
#' 
#' @return a data.frame corresponding to the C statistic, d.f., and P-value
#' 
#' @export
#'
fisherC <- function(dTable, add.claims = NULL, basis.set = NULL, direction = NULL, interactions = FALSE,
                    conserve = FALSE, conditional = FALSE, .progressBar = FALSE) {

  if(inherits(dTable, "list")) dTable <- as.psem(dTable)
  
  if(inherits(dTable, "psem")) dTable <- dSep(dTable, basis.set, direction, interactions, 
                                             conserve, conditional, .progressBar)

  if(length(dTable) == 0) {

    Cstat <- NA

    DF <- 0

    P <- NA

  } else {

    ps <- dTable$P.Value

    if(!is.null(add.claims)) {

      ps <- c(ps, add.claims)

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

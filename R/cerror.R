#' Correlated error operator
#'
#' Specifies correlated errors among predictors
#'
#' For use in \code{psem} to identify correlated sets of variables.
#'
#' @param e1 first variable involved in correlated error
#' @param e2 second variable involved in correlated error
#' 
#' @author Jon Lefcheck <LefcheckJ@@si.edu>, Jarrett Byrnes
#' @seealso \code{\link{cerror}}
#' 
#' @examples
#' # Generate example data
#' dat <- data.frame(x1 = runif(50),
#'   x2 = runif(50), y1 = runif(50),
#'     y2 = runif(50))
#'
#' # Create list of structural equations
#' sem <- psem(
#'   lm(y1 ~ x1 + x2, dat),
#'   lm(y2 ~ y1 + x1, dat)
#' )
#'
#' # Look at correlated error between x1 and x2
#' # (exogenous)
#' cerror(x1 %~~% x2, sem, dat)
#'
#' # Same as cor.test
#' with(dat, cor.test(x1, x2))
#'
#' # Look at correlatde error between x1 and y1
#' # (endogenous)
#' cerror(y1 %~~% x1, sem, dat)
#'
#' # Not the same as cor.test
#' # (accounts for influence of x1 and x2 on y1)
#' with(dat, cor.test(y1, x1))
#'
#' # Specify in psem
#' sem <- update(sem, x1 %~~% y1)
#'
#' coefs(sem)
#'
#' @export
#' 
`%~~%` <- function(e1, e2) {
  
  x <- paste(deparse(substitute(e1)), "~~", deparse(substitute(e2)))
  
  # x <- call(x)
  
  class(x) <- "formula.cerror"
  
  return(x)
  
}

#' Correlated errors
#'
#' Calculates partial correlations and partial significance tests.
#'
#' If the variables are exogenous, then the correlated error is the raw
#' bivariate correlation.
#'
#' If the variables are endogenous, then the correlated error is the partial
#' correlation, accounting for the influence of any predictors.
#'
#' The significance of the correlated error is conducted using \code{cor.test}
#' if the variables are exogenous. Otherwise, a t-statistic is constructed and
#' compared to a t-distribution with N - k - 2 degrees of freedom (where N is
#' the total number of replicates, and k is the total number of variables
#' informing the relationship) to derive a P-value.
#'
#' @param formula.  A formula specifying the two correlated variables using \code{\%~~\%}.
#' @param modelList A list of structural equations.
#' @param data A \code{data.frame} containing the data used in the list of equations.
#' 
#' @return Returns a \code{data.frame} containing the (partial) correlation and
#' associated significance test.
#' @author Jon Lefcheck <lefcheckj@@si.edu>
#' @seealso \code{\link{\%~~\%}}
#' 
#' @examples
#' # Generate example data
#' dat <- data.frame(x1 = runif(50),
#'   x2 = runif(50), y1 = runif(50),
#'     y2 = runif(50))
#'
#' # Create list of structural equations
#' sem <- psem(
#'   lm(y1 ~ x1 + x2, dat),
#'   lm(y2 ~ y1 + x1, dat)
#' )
#'
#' # Look at correlated error between x1 and x2
#' # (exogenous)
#' cerror(x1 %~~% x2, sem, dat)
#'
#' # Same as cor.test
#' with(dat, cor.test(x1, x2))
#'
#' # Look at correlatde error between x1 and y1
#' # (endogenous)
#' cerror(y1 %~~% x1, sem, dat)
#'
#' # Not the same as cor.test
#' # (accounts for influence of x1 and x2 on y1)
#' with(dat, cor.test(y1, x1))
#'
#' # Specify in psem
#' sem <- update(sem, x1 %~~% y1)
#'
#' coefs(sem)
#'
#' @export 
#' 
cerror <- function(formula., modelList, data = NULL) {
  
  ret <- partialCorr(formula., modelList, data)
  
  # ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)
  
  ret <- cbind.data.frame(ret, isSig(ret$P.Value))
  
  names(ret)[ncol(ret)] <- ""
  
  return(ret)
  
}

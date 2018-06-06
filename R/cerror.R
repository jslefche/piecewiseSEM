#' Correlated error operator
#'
#' Specifies correlated errors among predictors
#'
#' For use in \code{psem} to identify correlated sets of variables.
#' 
#' @usage e1 %~~% e2
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @aliases `~~`
#' @seealso \code{\link{cerror}}
#' @examples
#'
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
#' @export `%~~%`
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
#' @return Returns a \code{data.frame} containing the (partial) correlation and
#' associated significance test.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso \code{\link{\%~~\%}}
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
#' @export cerror
cerror <- function(formula., modelList, data = NULL) {

  tab <- partialCorr(formula., modelList, data)

  tab[, which(sapply(tab, is.numeric))] <- round(tab[, which(sapply(tab, is.numeric))], 4)

  return(tab)

}

#' Identify models with correlated errors and return modified versions
#' 
#' @keywords internal
#' 
getResidModels <- function(vars, modelList, data) {

  yvar <- sapply(listFormula(modelList), function(x) vars[[1]] %in% all.vars.merMod(x)[1])

  if(all(yvar == FALSE)) {

    vars <- rev(vars)

    yvar <- sapply(listFormula(modelList), function(x) vars[[1]] %in% all.vars.merMod(x)[1])

  }

  xvar <- sapply(listFormula(modelList), function(x) all(vars[[2]] %in% all.vars.merMod(x)[1]))

  if(all(yvar == FALSE) & all(xvar == FALSE)) {

    rdata <- data[, colnames(data) %in% vars]

    ymod <- as.numeric(data[, vars[[1]]])

    names(ymod) <- rownames(data)

    xmod <- as.numeric(data[, vars[[2]]])

    names(xmod) <- rownames(data)

  } else {

    if(all(xvar == FALSE)) {

      xvar <- sapply(listFormula(modelList), function(x) {

        f <- all.vars.merMod(x)

        any(f[1] == vars[[1]] & f[-1] %in% vars[[2]])

      } )

    }

    ymod <- modelList[[which(yvar)]]

    # if(length(all.vars.merMod) < 3) stop("Variables are part of a simple linear regression: partial residuals cannot be calculated!")

    termlabels.y <- which(grepl(paste(vars[[2]], collapse = ":"), all.vars_notrans(ymod)[-1]))

    if(length(termlabels.y) == 0) {

      vars[[2]] <- rev(vars[[2]])

      termlabels.y <- which(grepl(paste(vars[[2]], collapse = ":"), all.vars_notrans(ymod)[-1]))

    }

    if(length(termlabels.y) > 0) ymod <- update(ymod, drop.terms(terms(ymod), termlabels.y, keep.response = TRUE))

    if(all(xvar == FALSE)) {

      xmod <- as.numeric(data[, vars[[2]]])

      names(xmod) <- rownames(data)

    } else {

      xmod <- modelList[[which(xvar)]]

      newyvar <- all.vars_trans(xmod)[which(paste(vars[[2]], collapse = ":") == all.vars_notrans(xmod))]

      if(length(vars[[2]]) > 1) {

        splitxvar <- unlist(strsplit(newyvar, ":"))

        newdata <- data

        for(i in 1:length(splitxvar)) {

          newdata[, vars[[2]][i]] <- sapply(newdata[, vars[[2]][i]], function(x) eval(parse(text = gsub(vars[[2]][i], x, splitxvar[i]))))

        }

        newdata <- data.frame(newdata, apply(newdata[, vars[[2]]], 1, prod, na.rm = TRUE))

        data <- data.frame(data, newdata[, ncol(newdata)])

        names(data)[ncol(data)] <- paste(vars[[2]], collapse = "......")

        xmod <- update(xmod,
                       formula(paste(paste(vars[[2]], collapse = "......"), "~ ", paste(all.vars_trans(ymod)[-1], collapse = " + "))),
                       data = data)

      } else {

        if(length(termlabels.y) > 0) {

          f <- paste(newyvar, " ~ ", paste(all.vars_trans(ymod)[-1], collapse = " + "))

          xmod <- update(xmod, formula(f))

        }

      }

    }

  }

  list(ymod = ymod, xmod = xmod)

}

#' Correlated error operator
#'
#' Specifies correlated errors among predictors
#'
#' For use in \code{psem} to identify correlated sets of variables.
#'
#' @return Returns a formula of class \code{formula.cerror}.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
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
#' @export %~~%
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
#' @param formula.  A formula specifying the two correlated variables using
#' \code{%~~%}.
#' @param modelList A list of structural equations using \code{psem}.
#' @param data A \code{data.frame} containing the data used in the list of
#' equations.
#' @return Returns a \code{data.frame} containing the (partial) correlation and
#' associated significance test.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso \code{\link{cor.test}}, \code{\link{%~~%}}
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
#' @export cerror
cerror <- function(formula., modelList, data = NULL) {

  tab <- partialCorr(formula., modelList, data)

  tab[, which(sapply(tab, is.numeric))] <- round(tab[, which(sapply(tab, is.numeric))], 4)

  return(tab)

}

#' Identify models with correlated errors and return modified versions
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

    termlabels.y <- which(grepl(paste(vars[[2]], collapse = ":"), all.vars.notrans(ymod)[-1]))

    if(length(termlabels.y) == 0) {

      vars[[2]] <- rev(vars[[2]])

      termlabels.y <- which(grepl(paste(vars[[2]], collapse = ":"), all.vars.notrans(ymod)[-1]))

    }

    if(length(termlabels.y) > 0) ymod <- update(ymod, drop.terms(terms(ymod), termlabels.y, keep.response = TRUE))

    if(all(xvar == FALSE)) {

      xmod <- as.numeric(data[, vars[[2]]])

      names(xmod) <- rownames(data)

    } else {

      xmod <- modelList[[which(xvar)]]

      newyvar <- all.vars.trans(xmod)[which(paste(vars[[2]], collapse = ":") == all.vars.notrans(xmod))]

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
                       formula(paste(paste(vars[[2]], collapse = "......"), "~ ", paste(all.vars.trans(ymod)[-1], collapse = " + "))),
                       data = data)

      } else {

        if(length(termlabels.y) > 0) {

          f <- paste(newyvar, " ~ ", paste(all.vars.trans(ymod)[-1], collapse = " + "))

          xmod <- update(xmod, formula(f))

        }

      }

    }

  }

  list(ymod = ymod, xmod = xmod)

}

#' Computing partial effects
#'
#' Extracts partial residuals from a model or \code{psem} object for a given
#' \code{x} and \code{y}.
#'
#' This function computes the partial residuals of \code{y ~ x + Z} in a
#' two-step procedure to remove the variation explained by \code{Z}: (1) remove
#' \code{x} from the equation and model \code{y ~ Z}, and (2) replace \code{y}
#' with \code{x} and model \code{x ~ Z}.
#'
#' @param formula.  A formula where the \code{lhs} is the response and the
#' \code{rhs} is the predictor whose partial effect is desired.
#' @param modelList A list of structural equations.
#' @param data A \code{data.frame} used to fit the equations.
#' @return Returns a \code{data.frame} of residuals of \code{y ~ Z} called
#' \code{yresids}, of \code{x ~ Z} called \code{xresids}.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso \code{\link{cerror}}
#' @examples
#'
#' # Generate data
#' dat <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
#'
#' # Build model
#' model <- lm(y ~ x1 + x2, dat)
#'
#' # Compute partial residuals of y ~ x1
#' yresid <- resid(lm(y ~ x2, dat))
#'
#' xresid <- resid(lm(x1 ~ x2, dat))
#'
#' plot(yresid, xresid)
#'
#' # Use partialResid
#' presid <- partialResid(y ~ x1, model)
#'
#' plot(presid) # identical plot!
#'
#' @export partialResid
partialResid <- function(formula., modelList, data = NULL) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- removeData(modelList, formulas = 1)

  if(class(formula.) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(formula., "~~"))) else

      vars <- gsub(" ", "", unlist(strsplit(deparse(formula.), "~")))

  vars <- strsplit(vars, ":|\\*")

  if(!all(unlist(vars) %in% colnames(data)))

    stop("Variables not found in the model list. Ensure spelling is correct and remove all transformations!")

  residModList <- getResidModels(vars, modelList, data)

  if(all(class(residModList$ymod) == "numeric"))

    yresid <- data.frame(.id = names(residModList$ymod), yresid = residModList$ymod) else

      yresid <- data.frame(.id = rownames(getData.(residModList$ymod)), yresid = as.numeric(resid(residModList$ymod))) #resid.lme(ymod)

  if(all(class(residModList$xmod) == "numeric"))

    xresid <- data.frame(.id = names(residModList$xmod), xresid = residModList$xmod) else

      xresid <- data.frame(.id = rownames(getData.(residModList$xmod)), xresid = as.numeric(resid(residModList$xmod))) #resid.lme(xmod)

  rdata <- merge(yresid, xresid, by = ".id", all = TRUE)

  rdata <- rdata[order(as.numeric(as.character(rdata$.id))), -1]

  rownames(rdata) <- NULL

  return(rdata)

}

#' Calculate partial correlations from partial residuals
partialCorr <- function(formula., modelList, data = NULL) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- removeData(modelList, formulas = 1)

  rdata <- partialResid(formula., modelList, data)

  rcor <- cor(rdata[, 1], rdata[, 2], use = "complete.obs")

  if(class(formula.) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(formula., "~~"))) else

    vars <- gsub(" ", "", unlist(strsplit(deparse(formula.), "~")))

  vars <- strsplit(vars, ":|\\*")

  flag <- unlist(vars) %in% unlist(sapply(listFormula(modelList), function(x) all.vars.merMod(x)[1]))

  if(all(flag == FALSE)) {

    ctest <- cor.test(rdata[, 1], rdata[, 2])

    t. <- ctest$statistic

    N <- ctest$parameter

    P <- ctest$p.value

    } else {

      N <- nrow(rdata)

      residModList <- getResidModels(vars, modelList, data)

      k <- sum(sapply(residModList, function(x)

        if(all(class(x) == "numeric")) 0 else

          length(all.vars.merMod(formula(x))) - 2

      ) )

      k <- k[!duplicated(k)]

      k <- k[!k %in% vars]

      k <- length(k)

      N <- N - k - 2

      t. <- rcor * sqrt(N/(1 - rcor^2))

      P <- 1 - pt(t., N)

      }

  ret <- data.frame(
    Response = paste0("~~", vars[[1]]),
    Predictor = paste0("~~", paste(vars[[2]], collapse = ":")),
    Estimate = rcor,
    Std.Error = NA,
    DF = N,
    Crit.Value = t.,
    P.Value = P
  )

  return(ret)

}

#' Get residuals from innermost grouping of mixed models (replicate-level)
resid.lme <- function(model) {

  if(any(class(model) %in% c("lme", "glmmPQL"))) {

    Q <- length(summary(model)$modelStruct$reStruct)

    r <- resid(model, level = 0:Q)

    r <- r[, 1]

  } else r <- resid(model)

  return(r)

}

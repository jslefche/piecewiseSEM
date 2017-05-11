#' Warnings for old functions

partial.resid <- function(...) {

  warning("`partial.resid` has been replaced by `partialResid`", call. = FALSE)

}

sem.aic <- function(...) {

  warning("`sem.aic` has been replaced. Use `psem` instead of `list`, and then call `AIC` on that object", call. = FALSE)

}

sem.coefs <- function(...) {

  warning("`sem.coefs` has been replaced. Use `psem` instead of `list`, and then call `summary` or `coefs` on that object", call. = FALSE)

}

sem.fisher.c <- function(...) {

  warning("`sem.fisher.c` has been replaced. Use `psem` instead of `list`, and then call `summary` on that object", call. = FALSE)

}

sem.fit <- function(...) {

  warning("`sem.fit` has been replaced. Use `psem` instead of `list`, and then call `summary` on that object", call. = FALSE)

}

sem.model.fits <- function(...) {

  warning("`sem.model.fits` has been replaced by `rsquared`", call. = FALSE)

}

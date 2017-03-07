cerror <- function(.formula, modelList, data) {

  x <- pcor(.formula, modelList, data)

  p <- pcortest(x)

  return(x, p)

}

pcor <- function() {

}

pcortest <- function() {

}


`%~~%` <- function(e1, e2) {

  x <- paste(deparse(substitute(e1)), "~~", deparse(substitute(e2)))

  # x <- call(x)

  class(x) <- "formula.cerror"

  return(x)

}

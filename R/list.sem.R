#' Create a list of structural equations

#' @param ... arguments passed to methods

list.sem <- function(...) {


  # Evaluate correlated errors & store as new object

  # cerror()


  x <- list(...)

  evaluateClasses(x)

#   if(duplicated(sapply(listFormula(...), function(i) i[1])))
#
#     stop("Duplicated response variables detected. Collapse into single multiple regression and re-run.")
#
  class(x) <- "list.sem"

  return(print.attr(x))

}

# New operators for latent, composite variables


# `%~=%`

# `%~+%`

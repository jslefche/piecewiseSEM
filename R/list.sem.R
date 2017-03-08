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

  # remove quotes on cerrors
  print.attr(x)

  class(x) <- "sem"

  x

}

append.sem <- function(x, values) {

  x <- append(x, values)

  class(x[[length(x)]]) <- class(values)

  evaluateClasses(x)

  class(x) <- "sem"

  x

}



# New operators for latent, composite variables


# `%~=%`

# `%~+%`

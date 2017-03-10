#' Create a list of structural equations

#' @param ... arguments passed to methods

psem <- function(...) {


  # Evaluate correlated errors & store as new object

  # cerror()

  x <- list(...)

  evaluateClasses(x)

#   if(duplicated(sapply(listFormula(...), function(i) i[1])))
#
#     stop("Duplicated response variables detected. Collapse into single multiple regression and re-run.")
#

  # return error if duplicate responses
  if(any(duplicated(sapply(listFormula(x), function(y) all.vars.merMod(y)[1]))))

    stop("Duplicate responses detected in the model list. Collapse into single multiple regression!", call. = FALSE)

  # remove quotes on cerrors
  print.attr(x)

  class(x) <- "psem"

  x

}

# setMethod("list.sem", signature("sem"), function(...) list.sem(...))

update.psem <- function(x, values) {

  x <- append(x, values)

  class(x[[length(x)]]) <- class(values)

  evaluateClasses(x)

  class(x) <- "psem"

  x

}

# New operators for latent, composite variables

# `%~=%`

# `%~+%`

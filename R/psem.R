#' Create a list of structural equations

#' @param ... arguments passed to methods

psem <- function(...) {

  x <- list(...)

  evaluateClasses(x)

  formulaList <- listFormula(x)

  formulaList <- formulaList[!sapply(x, function(y) any(class(y) %in% c("formula", "formula.cerror")))]

  if(any(duplicated(sapply(formulaList, function(y) all.vars.merMod(y)[1]))))

    stop("Duplicate responses detected in the model list. Collapse into single multiple regression!", call. = FALSE)

  # remove quotes on cerrors
  # print.attr(x)

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

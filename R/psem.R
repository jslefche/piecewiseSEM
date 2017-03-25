#' Create a list of structural equations

#' @param ... arguments passed to methods

psem <- function(...) {

  x <- list(...)

  idx <- which(sapply(x, function(y) all(class(y) %in% c("matrix", "data.frame"))))

  x <- x[c((1:length(x))[!1:length(x) %in% idx], idx)]

  names(x)[length(x)] <- "data"

  evaluateClasses(x)

  formulaList <- listFormula(x)

  formulaList <- formulaList[!sapply(x, function(y) any(class(y) %in% c("matrix", "data.frame", "formula", "formula.cerror")))]

  if(any(duplicated(sapply(formulaList, function(y) all.vars.merMod(y)[1]))))

    stop("Duplicate responses detected in the model list. Collapse into single multiple regression!", call. = FALSE)

  # remove quotes on cerrors
  # print.attr(x)

  class(x) <- "psem"

  x

}

# setMethod("list.sem", signature("sem"), function(...) list.sem(...))

update.psem <- function(x, ...) {

  l <- list(...)

  for(i in 1:length(l)) x <- append(x, l[[i]])

  evaluateClasses(x)

  class(x) <- "psem"

  x

}

#' A list of supported model classes
model.classes <- c(
  "data.frame",
  "formula", "formula.cerror",
  "lm", "glm", "gls",
  "lme", "glmmPQL",
  "lmerMod", "merModLmerTest", "glmerMod"
  )

#' Evaluate model classes and stop if unsupported model class
evaluateClasses <- function(modelList) {

  classes <- unlist(sapply(modelList, class))

  classes <- classes[!duplicated(classes)]

  if(!all(classes %in% model.classes))

    stop(
      paste0(
        "Unsupported model class in model list: ",
        paste0(classes[!classes %in% model.classes], collapse = ", "),
        ". See 'help(piecewiseSEM)' for more details.")
    )

}

# New operators for latent, composite variables

# `%~=%`

# `%~+%`

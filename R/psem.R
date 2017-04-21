#' Create a list of structural equations
#'
#' @param ... arguments passed to methods

psem <- function(...) {

  x <- list(...)

  idx <- which(sapply(x, function(y) all(class(y) %in% c("matrix", "data.frame", "SpatialPointsDataFrame"))))

  if(length(idx) > 0) {

    x <- x[c((1:length(x))[!1:length(x) %in% idx], idx)]

    names(x)[length(x)] <- "data"

  } else {

    x$data <- getData.(x)

  }

  if(any(is.na(names(x)))) {

    idx. <- which(is.na(names(x)))

    names(x)[idx.] <- idx.

    }

  evaluateClasses(x)

  formulaList <- listFormula(x)

  formulaList <- formulaList[!sapply(x, function(y) any(class(y) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "formula", "formula.cerror")))]

  if(any(duplicated(sapply(formulaList, function(y) all.vars.merMod(y)[1]))))

    stop("Duplicate responses detected in the model list. Collapse into single multiple regression!", call. = FALSE)

  class(x) <- "psem"

  x

}

#' Update psem model object with additional values
update.psem <- function(x, ...) {

  l <- list(...)

  for(i in l) {

    if(all(class(i) %in% c("matrix", "data.frame"))) {

      idx <- which(names(x) == "data")

      if(length(idx) == 0) x$data = i else

        x[[idx]] <- i

      x <- lapply(x, function(j) {

        if(!any(class(j) %in% c("matrix", "data.frame", "formula", "formula.cerror")))

          update(j, data = i) else j

        } )

      } else if(all(class(i) %in% c("formula"))) {

        if(length(all.vars.merMod(i)) == 1) {

          x[[length(x) + 1]] <- i

        } else {

            resp <- sapply(x, function(y) if(!any(class(y) %in% c("matrix", "data.frame")))

              all.vars.merMod(y)[1] else "")

            idx <- which(resp == all.vars.merMod(i)[1])

            x[[idx]] <- update(x[[idx]], i)

            }

        } else {

          x[[length(x) + 1]] <- i

        }

  }

  evaluateClasses(x)

  class(x) <- "psem"

  x

}


#' Evaluate model classes and stop if unsupported model class
evaluateClasses <- function(modelList) {

  classes <- unlist(sapply(modelList, class))

  classes <- classes[!duplicated(classes)]

  model.classes <- c(
    "data.frame", "SpatialPointsDataFrame",
    "formula", "formula.cerror",
    "lm", "glm", "gls",
    "lme", "glmmPQL",
    "lmerMod", "merModLmerTest", "glmerMod",
    "sarlm",
    "pgls"
  )

  if(!all(classes %in% model.classes))

    stop(
      paste0(
        "Unsupported model class in model list: ",
        paste0(classes[!classes %in% model.classes], collapse = ", "),
        ". See 'help(piecewiseSEM)' for more details."),
      call. = FALSE
    )

}

# New operators for latent, composite variables

# `%~=%`

# `%~+%`

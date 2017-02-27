#' A list of supported model classes
model.classes <- c("noquote", "formula.cerror", "lm", "glm", "gls", "lme", "merMod", "glmerMod")

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

#' Get list of formula from a `sem` object
listFormula <- function(modelList) {

  # Get formulae
  x = lapply(modelList, function(i) if(any(class(i) == "formula.cerror")) i else formula(i) )

  # Strip transformations


  return(x)

}

#' Expand polynomials
expandPoly <- function(formulaList) {

  }

#' Strip transformations & offsets
stripTrans <- function(x) {

  }

#' Do not print attributes with custom functions
print.attr <- function(x) {

  attributes(x) <- NULL

  noquote(x)

}

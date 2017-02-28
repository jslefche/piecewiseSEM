#' A list of supported model classes
model.classes <- c("noquote", "formula.cerror", "lm", "glm", "gls", "lme", "lmerMod", "glmerMod")

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
#' If remove = TRUE, take out non-evaluated formula
listFormula <- function(modelList, remove = FALSE) {

  fList <- lapply(modelList, function(i) if(any(class(i) == "formula.cerror")) i else formula(i) )

  if(remove == TRUE) {

    l <- sapply(modelList, function(i) any(class(i) %in% c("formula.cerror")))

    fList <- fList[!l]

  }

  return(fList)

}

#' Remove random effects from all.vars
all.vars.merMod <- function(.formula) {

  n <- rownames(attr(terms(.formula), "factors"))

  if(any(grepl("\\|", n))) {

    idn <- which(grepl("\\|", n))

    f <- all.vars(.formula)[-idn]

    return(f)

  } else all.vars(.formula)

}

#' Expand polynomials
expandPoly <- function(formulaList) {

  }


#' Do not print attributes with custom functions
print.attr <- function(x) {

  attributes(x) <- NULL

  noquote(x)

}

#' Assess significance
isSig <- function(p) {

  ifelse(p > 0.01 & p < 0.05, "*",
       ifelse(p > 0.001 & p <= 0.01, "**",
              ifelse(p <= 0.001, "***", "")))

}

#' Get random effects from merMod
onlyBars <- function(.formula) {

  paste(

    sapply(lme4::findbars(.formula), function(x)

      paste0("(", deparse(x), ")")

    ),

    collapse = " + ")

}

#' Get list of formula from a `sem` object
#' If remove = TRUE, take out non-evaluated formula
listFormula <- function(modelList, remove = FALSE) {

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "formula")))]

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  fList <- lapply(modelList, function(i) if(any(class(i) %in% c("formula.cerror"))) i else formula(i) )

  if(remove == TRUE) {

    l <- sapply(modelList, function(i) any(class(i) %in% c("formula", "formula.cerror")))

    fList <- fList[!l]

  }

  return(fList)

}

#' Remove random effects from all.vars
all.vars.merMod <- function(.formula) {

  if(class(.formula) == "formula.cerror")

    gsub(" " , "", unlist(strsplit(.formula, "~~"))) else {

    n <- rownames(attr(terms(.formula), "factors"))

    if(any(grepl("\\|", n))) {

      f <- lme4::nobars(.formula)

      all.vars(f)

    } else all.vars(.formula)

    }

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

#' Recompute p-values using Kenward-Rogers approximation
KRp <- function(model, vars, intercepts = FALSE) {

  ret <- sapply(vars, function(x) {

    reducMod <- update(model, as.formula(paste("~ . -", x)))

    kr <- suppressWarnings(pbkrtest::KRmodcomp(model, reducMod))

    d <- round(kr$stats$ddf, 2)

    p <- kr$stats$p.valueU

    c(d, p)

  } )

  if(intercepts == TRUE) {

    reducMod <- update(model, as.formula(paste("~ 0 + .")))

    kr <- suppressWarnings(pbkrtest::KRmodcomp(model, reducMod))

    d <- kr$stats$ddf

    p <- kr$stats$p.valueU

    cbind(`(Intercept)` = c(d, p), ret)

  } else ret

}

#' Get data from model object
getData <- function(model) {

  if(any(class(model) %in% c("lm")))

    data <- model$model

  if(any(class(model) %in% c("glm", "glmmPQL")))

    data <- model$data

  if(any(class(model) %in% c("gls", "lme")))

    data <- nlme::getData(model)

  if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))

    data <- model@frame

  return(data)

}

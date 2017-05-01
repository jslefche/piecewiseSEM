#' Remove random effects from all.vars
all.vars.merMod <- function(formula.) {

  if(!any(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)

  if(class(formula.) == "formula.cerror")

    gsub(" " , "", unlist(strsplit(formula., "~~"))) else {

      n <- rownames(attr(terms(formula.), "factors"))

      if(any(grepl("\\|", n))) {

        f <- lme4::nobars(formula.)

        all.vars(f)

      } else all.vars(formula.)

    }

}

#' Get random effects from lme
findbars.lme <- function(model) {

  rand <- model$call$random

  rand <- gsub(".*\\|(.*)", "\\1", as.character(rand)[2])

  strsplit(gsub(" " , "", rand), "\\/")[[1]]

}

#' Get data from model list
getData. <- function(modelList) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  modelList <- removeData(modelList, formulas = 1)

  data.list <- lapply(modelList, function(model) {

    if(any(class(model) %in% c("lm", "sarlm")))

      data <- eval(model$call$data)

    if(any(class(model) %in% c("glm", "glmmPQL", "pgls")))

      data <- model$data

    if(any(class(model) %in% c("gls", "lme")))

      data <- nlme::getData(model)

    if(all(class(model) %in% c("pgls")))

      data <- model$data$data

    if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))

      data <- model@frame

    return(data)

  } )

  n <- max(sapply(data.list, nrow))

  data <- data.list[[which(sapply(data.list, function(x) nrow(x) == n))[1]]]

  for(i in 1:length(data.list)) {

    data <- merge(data, data.list[[i]], all = TRUE)

  }

  data <- data[, !duplicated(colnames(data), fromLast = TRUE)]

  return(data)

}

#' Get random effects variance-covariance from lme
getVarCov. <- function(model) {

  vc <- try(getVarCov(model), silent = TRUE)

  if(any(class(vc) == "try-error")) {

    vc <- VarCorr(model)

    v <- suppressWarnings(as.numeric(vc[, 1]))

    names(v) <- gsub(" =", "", rownames(vc))

    vm <- as.list(na.omit(v[-length(v)]))

    vl <- lapply(1:length(vm), function(i) matrix(vm[[i]], dimnames = list(names(vm)[i], names(vm)[i])))

    names(vl) <- names(which(is.na(v)))

    vl

  } else list(vc)

}

#' Assess significance
isSig <- function(p) {

  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))

}

#' Recompute P-values using Kenward-Rogers approximation
KRp <- function(model, vars, data, intercepts = FALSE) {

  ret <- sapply(vars, function(x) {

    reducMod <- update(model, as.formula(paste(". ~ . -", x)), data = data)

    kr <- suppressWarnings(pbkrtest::KRmodcomp(model, reducMod))

    d <- round(kr$stats$ddf, 2)

    p <- kr$stats$p.valueU

    c(d, p)

  } )

  if(intercepts == TRUE) {

    reducMod <- update(model, as.formula(paste("~ 0 + .")), data = data)

    kr <- suppressWarnings(pbkrtest::KRmodcomp(model, reducMod))

    d <- kr$stats$ddf

    p <- kr$stats$p.valueU

    cbind(`(Intercept)` = c(d, p), ret)

  } else ret

}

#' Get list of formula from a `sem` object
#' If remove = TRUE, take out non-evaluated formula
listFormula <- function(modelList, remove = FALSE) {

  modelList <- removeData(modelList, formulas = 0)

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  fList <- lapply(modelList, function(i) if(any(class(i) %in% c("formula.cerror"))) i else formula(i) )

  if(remove == TRUE) {

    l <- sapply(modelList, function(i) any(class(i) %in% c("formula", "formula.cerror")))

    fList <- fList[!l]

  }

  return(fList)

}

#' Get number of observations from a model
nobs. <- function(model) if(all(class(model) == "sarlm")) length(fitted(model)) else nobs(model)

#' Get random effects from merMod
onlyBars <- function(formula.) {

  paste(

    sapply(lme4::findbars(formula.), function(x)

      paste0("(", deparse(x), ")")

    ),

    collapse = " + ")

}

#' Do not print attributes with custom functions
print.attr <- function(x) {

  attributes(x) <- NULL

  noquote(x)

}

#' Remove data from the model list
removeData <- function(modelList, formulas = 0) {

  remove <- c("character", "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")

  if(formulas == 1) remove <- c(remove, "formula", "formula.cerror")

  if(formulas == 2) remove <- c(remove, "formula")

  modelList[!sapply(modelList, function(x) any(class(x) %in% remove))]

}

#' Remove random effects from all.vars
all.vars.merMod <- function(formula.) {

  if(!any(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)

  if(class(formula.) == "formula.cerror")

    gsub(" " , "", unlist(strsplit(formula., "~~"))) else {

      n <- rownames(attr(terms(formula.), "factors"))

      if(any(grepl("\\|", n)))

        all.vars(lme4::nobars(formula.)) else

          all.vars(formula.)

    }

}

#' Get vector of untransformed variables
all.vars.notrans <- function(formula.) {

  if(class(formula.) == "formula") {

    if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)

      all.vars(formula.)

    } else

      unlist(strsplit(formula., " ~~ "))

}




#' Get vector of transformed variables
all.vars.trans <- function(formula.) {

  if(class(formula.) == "formula") {

    if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)

    c(rownames(attr(terms(formula.), "factors"))[1], labels(terms(formula.)))

    } else unlist(strsplit(formula., " ~~ "))

}

#' Get random effects from lme
findbars.lme <- function(model) {

  rand <- model$call$random

  sapply(rand, function(i) {

    i = gsub(".*\\|(.*)", "\\1", as.character(i)[2])

    strsplit(gsub(" " , "", i), "\\/")[[1]]

  } )

}

#' Get data from model list
getData. <- function(modelList) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  modelList <- removeData(modelList, formulas = 1)

  data.list <- lapply(modelList, function(model) {

    if(any(class(model) %in% c("lm", "negbin", "sarlm")))

      data <- eval(model$call$data) else

    if(any(class(model) %in% c("glm", "glmmPQL")))

      data <- model$data else

    if(any(class(model) %in% c("gls", "lme")))

      data <- nlme::getData(model) else

    if(all(class(model) %in% c("pgls")))

      data <- model$data$data else

    if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))

      data <- model@frame else

        data <- data.frame()

    return(data)

  } )

  data <- Reduce(function(x, y) merge(x, y, all.x = FALSE), data.list)

  data <- data[, !duplicated(colnames(data), fromLast = TRUE)]

  colnames(data) <- gsub(".*\\((.*)\\).*", "\\1", colnames(data))

  data <- as.data.frame(data)

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

  if(grepl("\\*", deparse(formula(model))) & !all(grepl("\\*", vars))) {

    f <- all.vars.trans(formula(model))

    model <- update(model, as.formula(paste(f[1], " ~ ", paste(f[-1], collapse = " + "), " + ", paste(onlyBars(formula(model)), collapse = " + "))))

  }

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

    # cbind(`(Intercept)` = c(NA, NA), ret)

  } else ret

}

#' Get list of formula from a `sem` object
#' If remove = TRUE, take out non-evaluated formula
listFormula <- function(modelList, remove = FALSE) {

  modelList <- removeData(modelList, formulas = 0)

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  fList <- lapply(modelList, function(i) if(any(class(i) %in% c("formula.cerror"))) i else formula(i) )

  fList <- lapply(fList, lme4::nobars)

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
#' formulas = 0, keep everything
#' formulas = 1, remove all formulas including correlated errors
#' formulas = 2, remove only formula but keep correlated errors
removeData <- function(modelList, formulas = 0) {

  remove <- c("character", "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")

  if(formulas == 1) remove <- c(remove, "formula", "formula.cerror")

  if(formulas == 2) remove <- c(remove, "formula")

  modelList[!sapply(modelList, function(x) any(class(x) %in% remove))]

}

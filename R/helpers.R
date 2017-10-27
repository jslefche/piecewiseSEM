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

  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)

  if(class(formula.) == "formula") {

    if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)

    formula. <- all.vars.trans(formula.)

    if(any(grepl(":", formula.))) {

      idx <- which(grepl(":", formula.))

      for(i in idx) formula.[i] <- paste(sapply(strsplit(formula.[i], ":"), stripTransformations), collapse = ":")

      for(j in (1:length(formula.))[-idx]) formula.[j] <- stripTransformations(formula.[j])

    } else {

      formula. <- sapply(formula., stripTransformations)

    }

  } else formula. <- unlist(strsplit(formula., " ~~ "))

  return(formula.)

}

#' Get vector of transformed variables
all.vars.trans <- function(formula.) {

  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)

  if(class(formula.) == "formula") {

    if(formula.[[3]] == 1) deparse(formula.[[2]]) else {

      if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)

      c(rownames(attr(terms(formula.), "factors"))[1], labels(terms(formula.)))

      }

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

    if(any(class(model) %in% c("phylolm", "phyloglm"))) stop("Please provide data in `psem`")

    if(any(class(model) %in% c("lm", "negbin", "sarlm")))

      data <- eval(model$call$data) else

    if(any(class(model) %in% c("glm", "glmmPQL")))

      data <- model$data else

    if(any(class(model) %in% c("gls", "lme")))

      data <- nlme::getData(model) else

    if(all(class(model) %in% c("pgls"))) {

      data <- model$data

      }

    if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))

      data <- model@frame else

        data <- data.frame()

    return(data)

  } )

  if(all(sapply(data.list, class) == "comparative.data"))

    data <- data.list[[1]] else {

      data.list <- data.list[order(sapply(data.list, nrow))]

    if(length(data.list) > 1) {

      match.by <- unlist(sapply(data.list, names))

      match.by <- match.by[!duplicated(match.by)]

      data.list <- Map(function(x, i) setNames(x, ifelse(names(x) %in% match.by, names(x), sprintf('%s.%d', names(x), i))), data.list, seq_along(data.list))

      data <- Reduce(function(...) merge(..., all=T), data.list)

      } else data <- data.list[[1]]

    data <- data[, !duplicated(colnames(data), fromLast = TRUE)]

    colnames(data) <- gsub(".*\\((.*)\\).*", "\\1", colnames(data))

    data <- as.data.frame(data)

    }

  rownames(data) <- 1:nrow(data)

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
listFormula <- function(modelList, formulas = 0) {

  modelList <- removeData(modelList, formulas)

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  fList <- lapply(modelList, function(i) if(any(class(i) %in% c("formula.cerror"))) i else formula(i) )

  fList <- lapply(fList, lme4::nobars)

  return(fList)

}

#' Get number of observations from a model
nobs. <- function(object, ...) if(all(class(object) == "sarlm")) length(fitted(object)) else nobs(object, ...)

#' Get random effects from merMod
onlyBars <- function(formula., slopes = TRUE) {

  f <- lme4::findbars(formula.)

  if(slopes == TRUE) paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ") else {

    f <- f[sapply(f, function(x) grepl("1\\||1 \\|", deparse(x)))]

    paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ")

  }

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
#' formulas = 3, remove correlated errors but keep formula
removeData <- function(modelList, formulas = 0) {

  remove <- c("character", "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")

  if(formulas == 1) remove <- c(remove, "formula", "formula.cerror")

  if(formulas == 2) remove <- c(remove, "formula")

  if(formulas == 3) remove <- c(remove, "formula.cerror")

  modelList[!sapply(modelList, function(x) any(class(x) %in% remove))]

}

#' Strip transformations
stripTransformations <- function(x) {

  x <- gsub(".*\\((.*)\\).*", "\\1", x)

  gsub(" ", "", gsub("(.*)\\+.*", "\\1", x))

}

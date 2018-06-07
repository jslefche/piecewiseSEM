#' Remove random effects from all.vars
#' 
#' @keywords internal
#' 
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
#' 
#' @keywords internal
#' 
all.vars_notrans <- function(formula.) {

  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)

  if(class(formula.) == "formula") {

    if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)

    formula. <- all.vars_trans(formula.)

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
#' 
#' @keywords internal
#' 
all.vars_trans <- function(formula.) {

  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)

  if(class(formula.) == "formula") {

    if(formula.[[3]] == 1) deparse(formula.[[2]]) else {

      if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)

      c(rownames(attr(terms(formula.), "factors"))[1], labels(terms(formula.)))

      }

    } else unlist(strsplit(formula., " ~~ "))

}

#' Get random effects from lme
#' 
#' @keywords internal
#' 
findbars.lme <- function(model) {

  rand <- model$call$random

  sapply(rand, function(i) {

    i = gsub(".*\\|(.*)", "\\1", as.character(i)[2])

    strsplit(gsub(" " , "", i), "\\/")[[1]]

  } )

}

#' Get data from model list
#' 
#' @keywords internal
#' 
GetData <- function(modelList) {
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  modelList <- removeData(modelList, formulas = 1)
  
  data.list <- lapply(modelList, getSingleData)
  
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

#' Get data from one model
#' 
#' @keywords internal
#' 
getSingleData <- function(model) {

  dat <- data.frame()

  switch(class(model)[1],
         "phylolm" = {
           stop("Please provide `data =` argument to `psem`.", call. = FALSE)
         },

         "phyloglm" = {
           stop("Please provide `data =` argument to `psem`.", call. = FALSE)
         },

         "lm" ={
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },

         "negbin" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },

         "sarlm" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },

         "glm" = {
           dat <- model$data
         },

         "glmmPQL" = {
           dat <- model$data
         },

         "pgls" = {
           dat <- model$data
         },

         "lmerMod" = {
           dat <- model@frame
         },

         "glmerMod" = {
           dat <- model@frame
         },

         "merModLmerTest" = {
           dat <- model@frame
         },

         "gls" = {
           dat <-  nlme::getData(model)
         },

         "lme" = {
           dat <-  nlme::getData(model)
         }

  )

  dat

}

#' Obtain (observation-level) random effects from a generalized linear mixed model
#' 
#' RE = "all" all random effects are reported
#' RE = "RE" just group effects are reported
#' RE = "OLRE" just observation-level effects are reported
#' 
#' @keywords internal
#' 
GetOLRE <- function(sigma, model, X, data, RE = c("all", "RE", "OLRE")) {
  
  if(class(model) %in% c("lmerMod", "glmerMod")) {
    
    if(is.null(X)) X <- model.matrix(model)
    
    rand <- sapply(lme4::findbars(formula(model)), function(x) as.character(x)[3])
    
    rand <- rand[!duplicated(rand)] } else
      
      if(class(model) %in% c("lme", "glmmPQL")) {
        
        
      }
  
  idx <- sapply(sapply(strsplit(rand, "\\:"), function(x) gsub("\\(|\\)", "", x)), function(x) {
    
    length(unique(data[, x])) == nrow(data)
    
  } )
  
  sigma.names <- unlist(strsplit(names(sigma), "\\."))
  
  idx. <- sapply(sigma.names, function(x) !any(x %in% rand[idx]))
  
  if(RE == "RE") 
    
    sapply(sigma[idx.], function(i) {
      
      Z <- as.matrix(X[, rownames(i), drop = FALSE])
      
      sum(rowSums(Z %*% i) * Z) / nrow(X)
      
    } ) else if(RE == "OLRE") 
      
      sapply(sigma[idx], function(i) {
        
        Z <- as.matrix(X[, rownames(i), drop = FALSE])
        
        sum(rowSums(Z %*% i) * Z) / nrow(X)
        
      } ) else if(RE == "all")
        
        sapply(sigma, function(i) {
          
          Z <- as.matrix(X[, rownames(i), drop = FALSE])
          
          sum(rowSums(Z %*% i) * Z) / nrow(X)
          
        } )
  
}

#' Get random effects variance-covariance from lme
#' 
#' @keywords internal
#' 
GetVarCov <- function(model) {

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
#' 
#' @keywords internal
#' 
isSig <- function(p) {

  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))

}

#' Recompute P-values using Kenward-Rogers approximation
#' 
#' @keywrods internal
#' 
KRp <- function(model, vars, data, intercepts = FALSE) {

  if(grepl("\\*", deparse(formula(model))) & !all(grepl("\\*", vars))) {

    f <- all.vars_trans(formula(model))

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
#' 
#' If remove = TRUE, take out non-evaluated formula
#' 
#' @keywords internal
#' 
listFormula <- function(modelList, formulas = 0) {

  modelList <- removeData(modelList, formulas)

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  fList <- lapply(modelList, function(i) if(any(class(i) %in% c("formula.cerror"))) i else formula(i) )

  fList <- lapply(fList, lme4::nobars)

  return(fList)

}

#' Get number of observations from a model
#' 
#' @keywords internal
#' 
nObs <- function(object, ...) if(any(class(object) %in% c("phylolm", "phyloglm", "sarlm"))) length(fitted(object)) else nobs(object, ...)

#' Get random effects from merMod
#' 
#' @keywords internal
#' 
onlyBars <- function(formula., slopes = TRUE) {

  f <- lme4::findbars(formula.)

  if(slopes == TRUE) paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ") else {

    f <- f[sapply(f, function(x) grepl("1\\||1 \\|", deparse(x)))]

    paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ")

  }

}

#' Do not print attributes with custom functions
#' 
#' @keywords internal
#' 
#' @export
#' 
print.attr <- function(x, ...) {

  attributes(x) <- NULL

  noquote(x)

}

#' Remove data from the model list
#' 
#' formulas = 0, keep everything
#' formulas = 1, remove all formulas including correlated errors
#' formulas = 2, remove only formula but keep correlated errors
#' formulas = 3, remove correlated errors but keep formula
#' 
#' @keywords internal
#' 
removeData <- function(modelList, formulas = 0) {

  remove <- c("character", "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")

  if(formulas == 1) remove <- c(remove, "formula", "formula.cerror")

  if(formulas == 2) remove <- c(remove, "formula")

  if(formulas == 3) remove <- c(remove, "formula.cerror")

  modelList[!sapply(modelList, function(x) any(class(x) %in% remove))]

}

#' Strip transformations
#' 
#' @keywords internal
#' 
stripTransformations <- function(x) {

  x <- gsub(".*\\((.*)\\).*", "\\1", x)

  gsub(" ", "", gsub("(.*)\\+.*", "\\1", x))

}

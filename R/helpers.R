#' Remove random effects from all.vars
#' 
#' @keywords internal
#' 
all.vars.merMod <- function(formula.) {

  if(!any(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)

  if(inherits(formula., "formula.cerror"))

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
  
  if(inherits(formula., "formula")) {
    
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
  
  ret <- gsub("(,.*)", "", formula.)
  
  return(ret)
  
}

#' Get vector of transformed variables
#' 
#' @keywords internal
#' 
all.vars_trans <- function(formula., smoothed = FALSE) {
  
  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(inherits(formula., "formula")) {
    
    if(formula.[[3]] == 1) ret <- deparse(formula.[[2]]) else {
      
      if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)
      
      ret <- c(rownames(attr(terms(formula.), "factors"))[1], labels(terms(formula.)))
      
      if(smoothed == FALSE) ret <- gsub("(.*)\\,.*", "\\1", gsub("s\\((.*)\\).*", "\\1", ret)) 
      
      # else {
        
        # ret <- gsub("(s\\(.*),.*", "\\1", ret)
        # 
        # if(any(grepl("s\\(", ret))) ret <- sapply(ret, function(x) 
        #   ifelse(grepl("s\\(", x) & !grepl("\\)", x), paste0(x, ")"), x))
        
      # }
      
      # ret <- gsub("(,.*)", "", ret)
      
    }
    
    return(ret)
    
  } else unlist(strsplit(formula., " ~~ "))
  
}

#' Captures output table
#' 
#' @keywords internal
#' 
captureTable <- function(g, row.names = FALSE) {
  
  g1 <- capture.output(print(g, row.names = row.names))
  
  if(all(g1 == "data frame with 0 columns and 0 rows")) 
    
    g1 <- "No independence claims present. Tests of directed separation not possible."
  
  g1 <- paste0(g1, "\n")
  
  return(g1)
  
}

#' Bind data.frames of differing dimensions
#'
#' From: https://stackoverflow.com/a/31678079
#' 
#' @param ... data.frames to be bound, separated by commas
#' 
#'  @keywords internal
#'   
cbind_fill <- function(...) {
  
  nm <- list(...) 
  
  dfdetect <- grepl("data.frame|matrix", unlist(lapply(nm, function(cl) paste(class(cl), collapse = " ") )))
  
  vec <- data.frame(nm[!dfdetect])
  
  n <- max(sapply(nm[dfdetect], nrow)) 
  
  vec <- data.frame(lapply(vec, function(x) rep(x, n)))
  
  if (nrow(vec) > 0) nm <- c(nm[dfdetect], list(vec))
  
  nm <- lapply(nm, as.data.frame)
  
  do.call(cbind, lapply(nm, function (df1) 
    
    rbind(df1, as.data.frame(matrix(NA, ncol = ncol(df1), nrow = n-nrow(df1), dimnames = list(NULL, names(df1))))) )) 

  }

#' Transform variables based on model formula and store in new data frame
#' 
#' @keywords internal
#' 
dataTrans <- function(formula., data) {
  
  notrans <- all.vars.merMod(formula.)
  
  if(inherits(formula., "formula.cerror")) notrans <- gsub(".*\\((.*)\\)", "\\1", notrans)
  
  trans <- all.vars_trans(formula.)
  
  trans <- unlist(strsplit(trans, "\\:"))
  
  trans <- trans[!duplicated(trans)]
  
  if(any(grepl("scale\\(.*\\)", trans))) {
    # 
    # trans[which(grepl("scale(.*)", trans))] <- notrans[which(grepl("scale(.*)", trans))]
    # 
    warning("`scale` applied directly to variable. Use argument `standardize = TRUE` instead.", call. = FALSE)
    
  }
  
  if(any(!notrans %in% trans)) {
    
    for(k in 1:length(notrans)) {
      
      if(is.factor(data[, notrans[k]])) next else 
        
        if(grepl("scale(.*)", trans[k])) data[, notrans[k]] <- scale(data[, notrans[k]]) else
          
          data[, notrans[k]] <-
            
            sapply(data[, notrans[k]], function(x) eval(parse(text = gsub(notrans[k], x, trans[k]))))
      
    }
    
  }
  
  colnames(data) <- notrans
  
  return(data)
  
}

#' Get ANOVA results from `merMod`
#'
#' @keywords internal
#'
getAnova <- function(model, test.statistic = "F", test.type = "III") {

  if(inherits(model, "glmmTMB")) test.statistic = "Chisq"
  
  krp <- as.data.frame(car::Anova(model, test.statistic = test.statistic, type = test.type))

  ct <- summary(model)$coefficients
  
  colnames(ct)[2] <- "Std.Error"
  
  ret <- do.call(rbind, lapply(1:nrow(krp), function(i) {
    
    if(rownames(krp)[i] %in% rownames(ct)) {
      
      cbind.data.frame(
        ct[i, 1:2, drop = FALSE],
        DF = krp[i, 3], 
        Crit.Value = krp[i, 1], 
        P = krp[i, ncol(krp)],
        row.names = NULL
      )
      
    } else {
      
      data.frame(
        Estimate = NA,
        Std.Error = NA,
        DF = krp[i, 3], 
        Crit.Value = krp[i, 1], 
        P = krp[i, ncol(krp)]
      )
      
    }
    
  } ) )
  
  # ret <- cbind.data.frame(
  #   ct[, 1:2], 
  #   DF = krp[, 3], 
  #   Crit.Value = krp[, 1], 
  #   P = krp[, ncol(krp)]
  # )

  names(ret)[ncol(ret)] <- "Pr(>|t|)"
  
  rownames(ret) <- rownames(krp)

  return(ret)

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
  
  data.list <- lapply(modelList, GetSingleData)
  
  data.list <- data.list[!sapply(data.list, is.null)]
  
  data.list <- unname(data.list)
  
  if(all(sapply(data.list, class) == "comparative.data"))
    
    data <- data.list[[1]] else 
      
      data <- do.call(cbind_fill, data.list)
  
  data <- data[, !duplicated(colnames(data), fromLast = TRUE)]
      
  # colnames(data) <- gsub(".*\\((.*)\\).*", "\\1", colnames(data))

  data <- as.data.frame(data)
  
  rownames(data) <- 1:nrow(data)
  
  return(data)
  
}

#' Get data from one model
#' 
#' @keywords internal
#' 
GetSingleData <- function(model) {

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

         "Sarlm" = {
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
           dat <- lme4::getData(model) #model@frame
         },

         "glmerMod" = {
           dat <- lme4::getData(model) #model@frame
         },

         "lmerModLmerTest" = {
           dat <- lme4::getData(model) #model@frame
         },
         
         "glmmTMB" = {
           dat <- model$frame 
         },

         "gls" = {
           dat <- nlme::getData(model)
         },

         "lme" = {
           dat <- nlme::getData(model)
         },
         "gam" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         }

  )
 
  return(dat)

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
    
    rand <- rand[!duplicated(rand)] 
    
  } 
  
  # else if(class(model) %in% c("lme", "glmmPQL")) { }
  
  idx <- sapply(sapply(strsplit(rand, "\\:"), function(x) gsub("\\(|\\)", "", x)), function(x) {
    
    length(unique(data[, x])) == nrow(data)
    
  } )
  
  sigma.names <- unlist(names(sigma)) # unlist(strsplit(names(sigma), "\\."))
  
  idx. <- sapply(sigma.names, function(x) !any(x %in% rand[idx]))
  
  if(RE == "RE") 
    
    out <- sapply(sigma[idx.], function(i) {
      
      if(all(rownames(i) %in% colnames(X))) X. <- X else
        
        X. <- do.call(cbind, model.matrix(model, type = "randomListRaw")) 
      
      Z <- as.matrix(X.[, rownames(i), drop = FALSE])
      
      sum(rowSums(Z %*% i) * Z) / nrow(X.)
      
    } ) else if(RE == "OLRE") {
      
      if(all(idx == FALSE)) out <- 0 else {
        
        out <- sapply(sigma[idx], function(i) {
          
          Z <- as.matrix(X[, rownames(i), drop = FALSE])
          
          sum(rowSums(Z %*% i) * Z) / nrow(X)
          
        } ) } } else if(RE == "all")
          
          out <- sapply(sigma, function(i) {
            
            Z <- as.matrix(X[, rownames(i), drop = FALSE])
            
            sum(rowSums(Z %*% i) * Z) / nrow(X)
            
          } )
  
  if(length(out) == 0) out <- 0
  
  return(out)
  
}

#' Get random effects variance-covariance from lme
#' 
#' @keywords internal
#' 
GetVarCov <- function(model) {

  vc <- try(getVarCov(model), silent = TRUE)

  if(any(class(vc) == "try-error")) {

    vc <- nlme::VarCorr(model)

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
#' @keywords internal
#' 
# KRp <- function(model, vars, data, intercepts = FALSE) {
# 
#   # if(any(grepl("\\*", all.vars_notrans(formula(model)))) & !all(grepl("\\*", vars))) {
# 
#     f <- all.vars_trans(formula(model))
# 
#     model <- update(model, as.formula(paste(f[1], " ~ ", paste(f[-1], collapse = " + "), " + ", paste(onlyBars(formula(model)), collapse = " + "))))
# 
#   # }
# 
#   out <- data.frame()
#   
#   for(x in vars) { #sapply(vars, function(x) {
# 
#     reduceModel <- update(model, as.formula(paste(". ~ . -", x)))
# 
#     if(nobs(model) != nobs(reduceModel)) stop("Different sample sizes for `KRmodcomp`. Remove all NAs and re-run")
#     
#     kr <- try(pbkrtest::KRmodcomp(model, reduceModel), silent = TRUE)
# 
#     if(class(kr) == "try-error") 
#       
#       stop("Cannot obtain P-values from `lmerMod` using `pbkrtest::KRmodcopm`. Consider fitting using `nlme::lme`") else {
#         
#         d <- round(kr$stats$ddf, 2)
#   
#         p <- kr$stats$p.valueU
#   
#         out <- rbind(out, data.frame(d, p))
#         
#       }
# 
#   } # )
# 
#   if(intercepts == TRUE) {
# 
#     reduceModelI <- update(model, as.formula(paste("~ . - 1")), data = data)
# 
#     krI <- try(pbkrtest::KRmodcomp(model, reduceModelI), silent = TRUE)
#     
#     if(class(krI) == "try-error") 
#       
#       stop("Cannot obtain P-values from `lmerMod` using `pbkrtest::KRmodcomp`. Consider re-fitting using `nlme::lme`")else {
#         
#         dI <- krI$stats$ddf
#     
#         pI <- krI$stats$p.valueU
#     
#         out <- rbind(data.frame(d = dI, p = pI), out)
#         
#       }
#       
#   }
#   
#   return(out)
# 
# }

#' Get list of formula from a `psem` object
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
nObs <- function(object, ...) if(any(class(object) %in% c("phylolm", "phyloglm", "Sarlm"))) length(fitted(object)) else nobs(object, ...)

#' Get random effects from merMod
#' 
#' @keywords internal
#' 
onlyBars <- function(formula., slopes = TRUE) {

  f <- lme4::findbars(formula.)

  if(slopes == TRUE) paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ") else {

    # paste(sapply(f, function(x) paste0("(1 ", gsub(".*(\\|.*)", "\\1", f), ")")), collapse = "+")
    
    f <- f[sapply(f, function(x) grepl("1\\||1 \\|", deparse(x)))]

    paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ")

  }

}

#' Do not print attributes with custom functions
#' 
#' @keywords internal
#' 
#' @method print attr
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

#' Get Response Name as a Character
#' 
#' @keywords internal
#' 
get_response <- function(mod) {

    mod <- removeData(mod)
  
    f <- lapply(mod, formula)

      r <- lapply(f, function(x) x[[2]])
  
      return(as.character(r))

}

#' Get Left-hand side of formulae
#' 
#' @keywords internal
#' 
getLHS <- function(formulaList){
  sapply(formulaList, function(x) as.character(x[[2]]))
}

#' Get Right-hand side of formulae
#' 
#' @keywords internal
#' 
getRHS <- function(formulaList){
  rhs <- sapply(formulaList, function(x) all.vars(x)[-1])
  unique(do.call(c, rhs))
}

#' Operator for non-overlap in sets
#' 
#' @keywords internal
#' 
"%not_in%" <- function(x, y) x[!x %in% y]


#' Get a sorted psem object in DAG order
#' 
#' @description Takes a [psem] object, pulls out the
#' DAG, and then sorts the psem object into the order
#' of the DAG (from exogenous to terminal endogenous
#' variable) for use by other functions. Note: removes
#' correlated errors.
#'
#' @param object A fit [psem] object
#' @param keepdata Defaults to TRUE. Should the
#' data with the psem be included in the returned
#' object?
#'
#' @return A new [psem] object, without the data.
#' @export
getSortedPsem <- function(object, keepdata = TRUE){
  #first, remove data
  dat <- object$data
  object <- removeData(object, formulas = 1)
  
  #Now, get formulae
  formulaList <- listFormula(object)
  lhs <- getLHS(formulaList)
  
  names(object)<- lhs
  
  #sort some dags so we do things in the right order  
  object_dag <- getDAG(formulaList)
  sorted_dag <- sortDag(object_dag, formulaList)
  lhs_sorted <- colnames(sorted_dag)
  lhs_sorted <- lhs_sorted[which(lhs_sorted %in% lhs)]
  
  #Sort the object
  object <- object[lhs_sorted]
  
  #should we include the data?
  if(keepdata) object$data <- dat
  
  #return
  return(object)
}

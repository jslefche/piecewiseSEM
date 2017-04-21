#' Operator: correlated errors

`%~~%` <- function(e1, e2) {

  x <- paste(deparse(substitute(e1)), "~~", deparse(substitute(e2)))

  # x <- call(x)

  class(x) <- "formula.cerror"

  return(x)

}

#' Calculating (partial) correlations
#'
#' @param formula. a formula
#' @param modelList a list of structural equations

cerror <- function(formula., modelList, data = NULL) {

  tab <- partialCorr(formula., modelList, data)

  tab[, which(sapply(tab, is.numeric))] <- round(tab[, which(sapply(tab, is.numeric))], 4)

  return(tab)

}

#' Extract partial residuals
partialResid <- function(formula., modelList, data = NULL) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- removeData(modelList, formulas = 1)

  if(class(formula.) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(formula., "~~"))) else

      vars <- gsub(" ", "", unlist(strsplit(deparse(formula.), "~")))

  vars <- strsplit(vars, ":|\\*")

  residModList <- getResidModels(vars, modelList, data)

  if(all(class(residModList$ymod) == "numeric"))

    yresid <- data.frame(.id = names(residModList$ymod), yresid = residModList$ymod) else

      yresid <- data.frame(.id = rownames(getData.(residModList$ymod)), yresid = as.numeric(resid(residModList$ymod))) #resid.lme(ymod)

  if(all(class(residModList$xmod) == "numeric"))

    xresid <- data.frame(.id = names(residModList$xmod), xresid = residModList$xmod) else

      xresid <- data.frame(.id = rownames(getData.(residModList$xmod)), xresid = as.numeric(resid(residModList$xmod))) #resid.lme(xmod)

  rdata <- merge(yresid, xresid, by = ".id", all = TRUE)[, -1]

  return(rdata)

}

#' Calculate partial correlations from partial residuals
partialCorr <- function(formula., modelList, data = NULL) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- removeData(modelList, formulas = 1)

  rdata <- partialResid(formula., modelList, data)

  rcor <- cor(rdata[, 1], rdata[, 2], use = "complete.obs")

  if(class(formula.) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(formula., "~~"))) else

    vars <- gsub(" ", "", unlist(strsplit(deparse(formula.), "~")))

  vars <- strsplit(vars, ":|\\*")

  flag <- unlist(vars) %in% unlist(sapply(listFormula(modelList), function(x) all.vars.merMod(x)[1]))

  if(all(flag == FALSE)) {

    ctest <- cor.test(rdata[, 1], rdata[, 2])

    t. <- ctest$statistic

    N <- ctest$parameter

    P <- ctest$p.value

    } else {

      N <- nrow(rdata)

      residModList <- getResidModels(vars, modelList, data)

      k <- sum(sapply(residModList, function(x)

        if(all(class(x) == "numeric")) 0 else

          length(all.vars.merMod(formula(x))) - 2

      ) )

      k <- k[!duplicated(k)]

      k <- k[!k %in% vars]

      k <- length(k)

      N <- N - k - 2

      t. <- rcor * sqrt(N/(1 - rcor^2))

      P <- 1 - pt(t., N)

      }

  ret <- data.frame(
    Response = paste0("~~", vars[[1]]),
    Predictor = paste0("~~", paste(vars[[2]], collapse = ":")),
    Estimate = rcor,
    Std.Error = NA,
    DF = N,
    Crit.Value = t.,
    P.Value = P
  )

  return(ret)

}

#' Get residuals from innermost grouping of mixed models (replicate-level)
resid.lme <- function(model) {

  if(any(class(model) %in% c("lme", "glmmPQL"))) {

    Q <- length(summary(model)$modelStruct$reStruct)

    r <- resid(model, level = 0:Q)

    r <- r[, 1]

  } else r <- resid(model)

  return(r)

}

#' Identify models with correlated errors and return modified versions
getResidModels <- function(vars, modelList, data) {

  yvar <- sapply(listFormula(modelList), function(x) vars[[1]] %in% all.vars.merMod(x)[1])

  if(all(yvar == FALSE)) {

    vars <- rev(vars)

    yvar <- sapply(listFormula(modelList), function(x) vars[[1]] %in% all.vars.merMod(x)[1])

  }

  xvar <- sapply(listFormula(modelList), function(x) all(vars[[2]] %in% all.vars.merMod(x)[1]))

  if(all(yvar == FALSE) & all(xvar == FALSE)) {

    if(!all(unlist(vars) %in% colnames(data)))

      stop("Variables not found in the model list. Ensure spelling is correct!") else {

        rdata <- data[, colnames(data) %in% vars]

        ymod <- data[, vars[[1]]]

        names(ymod) <- rownames(data)

        xmod <- data[, vars[[2]]]

        names(xmod) <- rownames(data)

      }

  } else {

    if(all(xvar == FALSE)) {

      xvar <- sapply(listFormula(modelList), function(x) {

        f <- all.vars.merMod(x)

        any(f[1] == vars[[1]] & f[-1] %in% vars[[2]])

      } )

      }

    ymod <- modelList[[which(yvar)]]

    if(length(vars[[2]]) > 1) {

      ymod <- update(ymod,
                     formula(paste(". ~ . -", paste(vars[[2]], collapse = ":"))))

      } else {

      ymod <- update(ymod, formula(paste(". ~ . -", vars[[2]])))

      }

    if(all(xvar == FALSE)) {

      xmod <- data[, vars[[2]]]

      names(xmod) <- rownames(data)

    } else {

      xmod <- modelList[[which(xvar)]]

      if(length(vars[[2]]) > 1) {

        data <- data.frame(data, apply(data[, vars[[2]]], 1, prod, na.rm = TRUE))

        names(data)[ncol(data)] <- paste(vars[[2]], collapse = "...")

        xmod <- update(xmod,
                       formula(paste(paste(vars[[2]], collapse = "..."), "~ . -", paste(vars[[2]], collapse = ":"))),
                       data = data)

      } else {

        xmod <- update(xmod, formula(paste(vars[[2]], " ~ . -", vars[[2]])))

      }

    }

  }

  list(ymod = ymod, xmod = xmod)

}

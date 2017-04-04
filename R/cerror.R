#' Operator: correlated errors

`%~~%` <- function(e1, e2) {

  x <- paste(deparse(substitute(e1)), "~~", deparse(substitute(e2)))

  # x <- call(x)

  class(x) <- "formula.cerror"

  return(x)

}

#' Calculating (partial) correlations

#' @param .formula a formula
#' @param modelList a list of structural equations

cerror <- function(.formula, modelList, data = NULL) {

  tab <- partialCorr(.formula, modelList, data)

  tab[, which(sapply(tab, is.numeric))] <- round(tab[, which(sapply(tab, is.numeric))], 4)

  return(tab)

}

#' Extract partial residuals
partialResid <- function(.formula, modelList, data = NULL) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList[[1]])

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "formula", "formula.cerror")))]

  if(class(.formula) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(.formula, "~~"))) else

      vars <- gsub(" ", "", unlist(strsplit(deparse(.formula), "~")))

  vars <- strsplit(vars, ":|\\*")

  yvar <- sapply(listFormula(modelList), function(x) vars[[1]] %in% all.vars.merMod(x)[1])

  xvar <- sapply(listFormula(modelList), function(x) all(vars[[2]] %in% all.vars.merMod(x)[1]))

  if(all(yvar == FALSE) & all(xvar == FALSE)) {

    if(!all(unlist(vars) %in% colnames(data)))

      stop("Variables not found in the model list. Ensure spelling is correct!") else

        rdata <- data[, colnames(data) %in% vars]

  } else {

    if(all(xvar == FALSE))

      xvar <- sapply(listFormula(modelList), function(x) {

        f <- all.vars.merMod(x)

        any(f[1] == vars[[1]] & f[-1] %in% vars[[2]])

      } )

    ymod <- modelList[[which(yvar)]]

    xmod <- modelList[[which(xvar)]]

    if(length(vars[[2]]) > 1) {

      ymod <- update(ymod,
                     formula(paste(". ~ . -", paste(vars[[2]], collapse = ":"))))

      data <- data.frame(data, apply(data[, vars[[2]]], 1, prod, na.rm = TRUE))

      names(data)[ncol(data)] <- paste(vars[[2]], collapse = "...")

      xmod <- update(xmod,
                     formula(paste(paste(vars[[2]], collapse = "..."), "~ . -", paste(vars[[2]], collapse = ":"))),
                     data = data)

      } else {

        ymod <- update(ymod, formula(paste(". ~ . -", vars[[2]])))

        xmod <- update(xmod, formula(paste(vars[[2]], " ~ . -", vars[[2]])))

        }

    yresid <- data.frame(.id = rownames(getData.(ymod)), yresid = resid(ymod)) #resid.lme(ymod)

    xresid <- data.frame(.id = rownames(getData.(xmod)), xresid = resid(xmod)) #resid.lme(xmod)

    rdata <- merge(yresid, xresid, by = ".id", all = TRUE)[, -1]

  }

  return(rdata)

}

#' Calculate partial correlations from partial residuals
partialCorr <- function(.formula, modelList, data = NULL) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList[[1]])

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "formula", "formula.cerror")))]

  rdata <- partialResid(.formula, modelList, data)

  rcor <- cor(rdata[, 1], rdata[, 2], use = "complete.obs")

  if(class(.formula) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(.formula, "~~"))) else

    vars <- gsub(" ", "", unlist(strsplit(deparse(.formula), "~")))

  vars <- strsplit(vars, ":|\\*")

  flag <- unlist(vars) %in% unlist(sapply(listFormula(modelList), function(x) all.vars.merMod(x)[1]))

  if(all(flag == FALSE)) {

    ctest <- cor.test(rdata[, 1], rdata[, 2])

    t. <- ctest$statistic

    N <- ctest$parameter

    P <- ctest$p.value

    } else {

      N <- nrow(rdata)

      k <- c(all.vars.merMod(formula(ymod)), all.vars.merMod(formula(xmod)))

      k <- k[!duplicated(k)]

      k <- k[!k %in% vars]

      k <- length(k)

      N <- N - k - 2

      t. <- rcor * sqrt(N/(1 - rcor^2))

      P <- pt(t., N)

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

#' Name vector so every entry is unique, preserving original order
name.vec <- function(v) {

  v <- sapply(unique(v), function(x) v[v == x] <- paste(v[v == x], 1:length(v[v == x]), sep = ".") )

  unlist(v)

}

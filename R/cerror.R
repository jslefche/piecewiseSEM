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

cerror <- function(.formula, modelList) {

  if(class(.formula) == "formula.cerror") {

    x <- pcor(.formula, modelList)

  }

  return(x)

}

pcor <- function(.formula, modelList) {

  if(class(.formula) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(.formula, "~~"))) else

      vars <- gsub(" ", "", unlist(strsplit(deparse(.formula), "~")))

  modelList2 <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("formula.cerror")))]

  yvar <- sapply(listFormula(modelList2), function(x) all.vars.merMod(x)[1] == vars[1])

  xvar <- sapply(listFormula(modelList2), function(x) all.vars.merMod(x)[1] == vars[2])

  if(all(yvar == FALSE) | all(xvar == FALSE))

    stop("Variables not found as responses in model list. Ensure spelling is correct!")

  ymod <- modelList[[which(yvar)]]

  if(!vars[2] %in% all.vars.merMod(formula(ymod))) {

    xmod <- modelList[[which(xvar)]]

    yresid <- resid.lme(ymod)

    yresid <- data.frame(.id = names(yresid), yresid = yresid)

    xresid <- resid.lme(xmod)

    xresid <- data.frame(.id = names(xresid), xresid = xresid)

    rdata <- merge(yresid, xresid, by = ".id", all = TRUE)[, -1]

  } else {


    # ymod <- update(

  }

  rcor <- cor(rdata[, 1], rdata[, 2], use = "complete.obs")

  N <- nrow(rdata)

  P <- 1 - pt((rcor * sqrt(N - 2))/(sqrt(1 - rcor^2)), (N - 2))

  ret <- data.frame(
    Response = paste0("~~", vars[1]),
    Predictor = paste0("~~", vars[2]),
    Estimate = rcor,
    Std.Error = NA,
    DF = N,
    Crit.Value = NA,
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

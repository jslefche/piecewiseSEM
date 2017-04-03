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

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList[[1]])

  if(class(.formula) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(.formula, "~~"))) else

      vars <- gsub(" ", "", unlist(strsplit(deparse(.formula), "~")))

  vars <- strsplit(vars, ":|\\*")

  yvar <- sapply(listFormula(modelList), function(x) all.vars.merMod(x)[1] %in% vars[[1]])

  xvar <- sapply(listFormula(modelList), function(x) any(all.vars.merMod(x)[-1] %in% vars[[2]]))

  if(all(yvar == FALSE) & all(xvar == FALSE)) {

    rdata <- data[, colnames(data) %in% vars]

  } else {

    if(which(yvar) != which(xvar))

      stop("Variables not found as responses in model list. Ensure spelling is correct!")

    ymod <- modelList[[which(yvar)]]

    if(length(vars[[2]]) > 1) {

      ymod <- update(ymod, formula(paste(". ~ . -", paste(vars[[2]], collapse = ":"))))

      data[, (ncol(data) + 1)] <- apply(data[, vars[[2]]], 1, prod, na.rm = TRUE)

      names(data)[ncol(data)] <- paste(vars[[2]], collapse = ".")

      xmod <- update(ymod, formula(paste(paste(vars[[2]], collapse = "."), "~ . -", paste(vars[[2]], collapse = ":"))))

      } else {

        ymod <- update(ymod, formula(paste(". ~ . -", vars[[2]])))

        xmod <- update(ymod, formula(paste(vars[[2]], " ~ . -", vars[[2]])))

        }

    yresid <- resid.lme(ymod)

    yresid <- data.frame(.id = names(yresid), yresid = yresid)

    xresid <- resid.lme(xmod)

    xresid <- data.frame(.id = names(xresid), xresid = xresid)

    rdata <- merge(yresid, xresid, by = ".id", all = TRUE)[, -1]

    }

  rcor <- cor(rdata[, 1], rdata[, 2], use = "complete.obs")

  if(all(yvar == FALSE) & all(xvar == FALSE)) {

    ctest <- cor.test(rdata[, 1], rdata[, 2])

    N <- ctest$parameter

    P <- ctest$p.value

    } else {

      N <- nrow(rdata)

      P <- 1 - pt((rcor * sqrt(N - 2))/(sqrt(1 - rcor^2)), (N - 2))

      }

  ret <- data.frame(
    Response = paste0("~~", vars[[1]]),
    Predictor = paste0("~~", paste(vars[[2]], collapse = ":")),
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

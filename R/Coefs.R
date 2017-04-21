#' Get (standardized) coefficients from list of structural equations
#'
#' @param modelList a list of structural equations

coefs <- function(modelList, data = NULL, intercepts = FALSE, standardize = TRUE) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- psem(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(class(data) %in% c("SpatialPointsDataFrame")) data <- data@data

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "formula")))]

  if(standardize == TRUE) ret <- stdCoefs(modelList, data, intercepts) else

    ret <- unstdCoefs(modelList, data, intercepts)

  ret <- cbind(ret, isSig(ret$P.Value))

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)

  names(ret)[length(ret)] = ""

  return(ret)

}

#' Get raw (understandardized) coefficients from model
unstdCoefs <- function(modelList, data = NULL, intercepts = FALSE) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "formula")))]

  ret <- do.call(rbind, lapply(modelList, function(i) {

    if(all(class(i) %in% c("formula.cerror")))

      ret <- cerror(i, modelList, data) else {

        if(all(class(i) %in% c("lm", "glm", "lmerMod", "glmerMod", "pgls"))) {

          ret <- as.data.frame(summary(i)$coefficients)

          if(all(class(i) %in% c("lm", "glm"))) ret <- cbind(ret[, 1:2], DF = summary(i)$df[2], ret[, 3:4])

          if(all(class(i) %in% c("glmerMod", "pgls"))) ret <- cbind(ret[, 1:2], DF = length(summary(i)$residuals), ret[, 3:4])

        }

        if(all(class(i) %in% c("sarlm"))) {

          ret <- as.data.frame(summary(i)$Coef)

          ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])

        }

        if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {

          krp <- KRp(i, all.vars.merMod(formula(i))[-1], intercepts = TRUE)

          ret <- as.data.frame(append(as.data.frame(ret), list(DF = krp[1,]), after = 2))

          ret[, "Pr(>|t|)"] <- krp[2, ]

          }

        if(all(class(i) %in% c("gls", "lme", "glmmPQL"))) {

          ret <- as.data.frame(summary(i)$tTable)

          if(ncol(ret) == 4 & class(i) %in% c("gls"))

            ret <- cbind(ret[, 1:2], DF = length(residuals(i)), ret[, 3:4])

          if(ncol(ret) == 4 & class(i) %in% c("lme", "glmmPQL"))

            ret <- cbind(ret[, 1:2], DF = summary(i)$fixDF$X[-1], ret[, 3:4])

          }

        ret <- data.frame(
          Response = all.vars.merMod(listFormula(list(i))[[1]])[1],
          Predictor = rownames(ret),
          ret
          )

        names(ret) <- c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value")

        }

    if(intercepts == FALSE) ret <- subset(ret, Predictor != "(Intercept)")

    return(ret)

    } ) )

  rownames(ret) <- NULL

  return(ret)

}

#' Calculate standardized regression coefficients
stdCoefs <- function(modelList, data = NULL, intercepts = FALSE) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "SpatialPointsDataFrame", "formula")))]

  ret <- unstdCoefs(modelList, data, intercepts)

  Bnew <- do.call(c, lapply(1:length(modelList), function(i) {

    j <- modelList[[i]]

    f <- listFormula(list(j))[[1]]

    if(all(all.vars.merMod(f) %in% colnames(data)))

      newdata <- data[, all.vars.merMod(f)] else

        newdata <- data[, all.vars.trans(f)]

    f <- all.vars.merMod(f)

    if(all(class(j) %in% c("formula.cerror"))) {

      Bnew <- subset(ret, Response == paste0("~~", f[1]) & Predictor == paste0("~~", f[2]))$Estimate

    } else {

      if(class(newdata) %in% c("SpatialPointsDataFrame")) newdata <- newdata@data

      newdata <- dataTrans(formula(j), newdata)

      B <- subset(ret, Response == f[1])$Estimate

      sd.x <- sapply(f[-1], function(x) sd(newdata[, x], na.rm = TRUE))

      if(length(B) != length(sd.x)) sd.x <- c(sd.x, sdInt(j, newdata))

      sd.y <- sdFam(f[1], j, newdata)

      if(intercepts == FALSE) Bnew <- B * (sd.x / sd.y) else

        Bnew <- c(0, B[-1] * (sd.x / sd.y))

    }

  } ) )

  cbind(ret, Std.Estimate = unname(Bnew))

}

#' Properly scale standard deviations depending on the error distribution
sdFam <- function(x, model, data) {

  .family = try(family(model), silent = TRUE)

  if(class(.family) == "try-error") y <- data[, x] else {

    if(.family$family != "gaussian")

      warning(
        paste0("Standardized coefficients for non-normal distributions must be interpreted differently. See help('stdCoefs')."),
        call. = FALSE
        )

    .link <- .family$link

    if(.link == "identity")

      y <- data[, x] else {

        if(any(class(model) %in% c("glmerMod")))

          linkfun <- model@resp$family$linkfun else

            linkfun <- model$family$linkfun

        y <- data[, x] * linkfun(mean(data[, x]))

      }

  }

  sd(y, na.rm = TRUE)

}

#' Transform variables based on model formula and store in new data frame
dataTrans <- function(formula., newdata) {

  notrans <- all.vars.notrans(formula.)

  trans <- all.vars.trans(formula.)

  if(any(grepl("scale\\(.*\\)", trans))) {

    trans[which(grepl("scale(.*)", trans))] <- notrans[which(grepl("scale(.*)", trans))]

    warning("`scale` applied to variable--use argument `standardize = TRUE` instead.", call. = FALSE)

  }

  if(any(notrans != trans)) {

    for(k in 1:length(notrans)) {

      newdata[, notrans[k]] <-

        sapply(newdata[, notrans[k]], function(x) eval(parse(text = gsub(notrans[k], x, trans[k]))))

    }

  }

  colnames(newdata) <- notrans

  return(newdata)

}

#' Calculate standard deviations for interactions
sdInt <- function(model, newdata) {

  v <- attr(terms(model), "term.labels")

  int <- v[grepl(":", v)]

  sapply(int, function(x) {

    sd(apply(newdata[, strsplit(x, ":")[[1]]], 1, prod, na.rm = TRUE))

  } )

}

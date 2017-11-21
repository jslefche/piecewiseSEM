#' Extract path coefficients
#'
#' Extracts (standardized) path coefficients from a \code{psem} object.
#'
#' P-values for models constructed using \code{\link{lme4::lme4}} are obtained
#' using the Kenward-Roger approximation of the denominator degrees of freedom
#' as implemented in the \code{\link{pbkrtest::pbkrtest}} package.
#'
#' Standardized coefficients are obtained by multiplying the parameter estimates
#' by the ratio of the standard deviation of x over the standard deviation of y.
#'
#' For binary response models (i.e., binomial responses), standardized coefficients
#' are obtained in one of two ways:\itemize{
#' \item{\code{latent} the latent theoretic approach,
#' which is the square-root of the variance of the predictions (on the linear or
#' 'link' scale) plus the distribution-specific variation (for logit links: pi^2/3, for
#' probit links: 1).}
#' \item{\code{observation} the observation error approach, where the square-root of the
#' variance of the predictions (on the linear scale) are divided by the correlation
#' between the observed and predicted (on th original or 'response' scale) values of y. }
#' }
#'
#' @param modelList A list of structural equations.
#' @param standardize Whether standardized coefficients should be reported.
#' Default is TRUE.
#' @param standardize.type The type of standardized for non-Gaussian responses.
#' Current options are \code{latent} or \code{observation}. Default is \code{observation}.
#' @param intercepts Whether intercepts should be included in the coefficients
#' table. Default is FALSE.
#' @return Returns a \code{data.frame} of coefficients, their standard errors,
#' degrees of freedom, and significance tests.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso \code{\link{KRmodcomp}}
#' @examples
#' @export coefs
coefs <- function(modelList, standardize = TRUE, standardize.type = "observation", intercepts = FALSE) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(class(modelList) == "psem") data <- modelList$data else data <- getData.(modelList)

  if(class(data) %in% c("SpatialPointsDataFrame")) data <- data@data

  if(class(data) %in% c("comparative.data")) data <- data$data

  if(standardize == TRUE) ret <- stdCoefs(modelList, data, standardize.type, intercepts) else

    ret <- unstdCoefs(modelList, data, intercepts)

  ret <- cbind(ret, isSig(ret$P.Value))

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)

  names(ret)[length(ret)] = ""

  return(ret)

}

#' Get raw (understandardized) coefficients from model
unstdCoefs <- function(modelList, data = NULL, intercepts = FALSE) {

  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- removeData(modelList, formulas = 2)

  ret <- do.call(rbind, lapply(modelList, function(i) {

    if(all(class(i) %in% c("formula.cerror")))

      ret <- cerror(i, modelList, data) else {

        if(all(class(i) %in% c("lm", "glm", "negbin", "lmerMod", "glmerMod", "merModLmerTest", "pgls", "phylolm", "phyloglm"))) {

          ret <- as.data.frame(summary(i)$coefficients)

          if(all(class(i) %in% c("lm", "glm", "negbin"))) ret <- cbind(ret[, 1:2], DF = summary(i)$df[2], ret[, 3:4])

          if(all(class(i) %in% c("glmerMod", "pgls"))) ret <- cbind(ret[, 1:2], DF = length(summary(i)$residuals), ret[, 3:4])

          if(all(class(i) %in% c("phylolm", "phyloglm"))) ret <- cbind(ret[, 1:2], DF = bNewMod$n, ret[, c(3, 6)])

          if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {

            krp <- KRp(i, all.vars.trans(formula(i))[-1], data, intercepts = TRUE)

            ret <- as.data.frame(append(as.data.frame(ret), list(DF = krp[1,]), after = 2))

            ret[, "Pr(>|t|)"] <- krp[2, ]

          }

        }

        if(all(class(i) %in% c("sarlm"))) {

          ret <- as.data.frame(summary(i)$Coef)

          ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])

        }

        if(all(class(i) %in% c("gls", "lme", "glmmPQL"))) {

          ret <- as.data.frame(summary(i)$tTable)

          if(ncol(ret) == 4 & all(class(i) %in% c("gls")))

            ret <- cbind(ret[, 1:2], DF = length(residuals(i)), ret[, 3:4])

          }

        ret <- data.frame(
          Response = all.vars.trans(listFormula(list(i))[[1]])[1],
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
stdCoefs <- function(modelList, data = NULL, standardize.type = "observation", intercepts = FALSE) {

  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- removeData(modelList, formulas = 2)

  ret <- unstdCoefs(modelList, data, intercepts)

  Bnew <- do.call(rbind, lapply(1:length(modelList), function(i) {

    j <- modelList[[i]]

    f <- listFormula(list(j))[[1]]

    newdata <- data[, all.vars.merMod(f)]

    f.trans <- all.vars.trans(f)

    f.notrans <- all.vars.notrans(f)

    if(all(class(j) %in% c("formula.cerror"))) {

      Bnew <- subset(ret, Response == paste0("~~", f.trans[1]) & Predictor == paste0("~~", f.trans[2]))$Estimate

    } else {

      if(any(class(newdata) %in% c("SpatialPointsDataFrame"))) newdata <- newdata@data

      newdata. <- dataTrans(formula(j), newdata)

      B <- subset(ret, Response == f.trans[1])$Estimate

      B.se <- subset(ret, Response == f.trans[1])$Std.Error

      sd.x <- sapply(f.notrans[!grepl(":", f.notrans)][-1], function(x) sd(newdata.[, x], na.rm = TRUE))

      if(any(grepl(":", f.notrans))) sd.x <- c(sd.x, sdInt(j, newdata.))

      sd.y <- sdFam(f.notrans[1], j, newdata., standardize.type)

      if(intercepts == FALSE)

        data.frame(Std.Estimate = B * (sd.x / sd.y), Std.SE = B.se * (sd.x / sd.y)) else

          data.frame(Std.Estimate = c(0, B[-1] * (sd.x / sd.y)), Std.SE = c(0, B.se[-1] * (sd.x / sd.y)))

      }

  } ) )

  ret <- data.frame(ret, Bnew)

  rownames(ret) <- NULL

  return(ret)

}

#' Transform variables based on model formula and store in new data frame
dataTrans <- function(formula., newdata) {

  notrans <- all.vars.merMod(formula.)

  trans <- all.vars.trans(formula.)

  trans <- unlist(strsplit(trans, "\\:"))

  trans <- trans[!duplicated(trans)]

  if(any(grepl("scale\\(.*\\)", trans))) {

    trans[which(grepl("scale(.*)", trans))] <- notrans[which(grepl("scale(.*)", trans))]

    warning("`scale` applied directly to variable. Use argument `standardize = TRUE` instead.", call. = FALSE)

  }

  if(any(!notrans %in% trans)) {

    for(k in 1:length(notrans)) {

      newdata[, notrans[k]] <-

        sapply(newdata[, notrans[k]], function(x) eval(parse(text = gsub(notrans[k], x, trans[k]))))

    }

  }

  colnames(newdata) <- notrans

  return(newdata)

}

#' Properly scale standard deviations depending on the error distribution
sdFam <- function(x, model, newdata, standardize.type = "observation") {

  .family <- try(family(model), silent = TRUE)

  if(class(.family) == "try-error") .family <- try(model$family, silent = TRUE)

  if(class(.family) == "try-error") NA else {

    if(.family$family == "gaussian") sd.y <- sd(newdata[, x], na.rm = TRUE)

    if(.family$family == "binomial") sd.y <- sdGLM(model, standardize.type)

    # ... additional model types here

  }

  return(sd.y)

}

#' Compute standard deviation of response for GLMs
sdGLM <- function(model, standardize.type = "observation") {

  preds <- predict(model, type = "link")

  if(standardize.type == "observation") {

    y <- all.vars.notrans(model)[1]

    data <- getSingleData(model)

    R <- cor(data[, y], predict(model, type = "response"))

    sd.y <- sqrt(var(preds)) / R

  }

  if(standardize.type == "latent") {

    link. <- family(model)$link

    if(link. == "logit") sigmaE <- pi^2/3 else

      if(link. == "probit") sigmaE <- 1

    sd.y <- sqrt(var(preds) + sigmaE)

  }

  return(sd.y)

}

#' Calculate standard deviations for interactions
sdInt <- function(model, newdata) {

  v <- attr(terms(model), "term.labels")

  int <- v[grepl(":", v)]

  sapply(int, function(x) {

    x <- strsplit(x, ":")[[1]]

    x <- gsub("(.*) \\+.*", "\\1", gsub(".*\\((.*)\\)", "\\1", x))

    sd(apply(newdata[, x], 1, prod, na.rm = TRUE))

  } )

}

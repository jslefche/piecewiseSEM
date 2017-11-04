#' Extract path coefficients
#' 
#' Extracts (standardized) path coefficients from a \code{psem} object.
#' 
#' P-values for models constructed using \code{\link{lme4::lme4}} are obtained
#' using the Kenward-Roger approximation of the denominator degrees of freedom
#' as implemented in the \code{\link{pbkrtest::pbkrtest}} package.
#' 
#' @param modelList A list of structural equations.
#' @param intercepts Whether intercepts should be included in the coefficients
#' table. Default is FALSE.
#' @param standardize Whether standardized path coefficients should be
#' reported. Default is TRUE.  ~~Describe \code{standardize} here~~
#' @return Returns a \code{data.frame} of coefficients, their standard errors,
#' degrees of freedom, and significance tests.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso \code{\link{KRmodcomp}}
#' @examples
#' 
#' 
#' 
#' @export coefs
coefs <- function(modelList, intercepts = FALSE, standardize = TRUE) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(class(modelList) == "psem") data <- modelList$data else data <- getData.(modelList)

  if(class(data) %in% c("SpatialPointsDataFrame")) data <- data@data

  if(class(data) %in% c("comparative.data")) data <- data$data

  modelList <- removeData(modelList, formulas = 2)

  if(standardize == TRUE) ret <- stdCoefs(modelList, data, intercepts) else

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
stdCoefs <- function(modelList, data = NULL, intercepts = FALSE) {

  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)

  if(is.null(data) & class(modelList) == "psem") data <- modelList$data

  if(is.null(data)) data <- getData.(modelList)

  modelList <- removeData(modelList, formulas = 2)

  ret <- unstdCoefs(modelList, data, intercepts)

  Bnew <- do.call(c, lapply(1:length(modelList), function(i) {

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

      sd.x <- sapply(f.notrans[!grepl(":", f.notrans)][-1], function(x) sd(newdata.[, x], na.rm = TRUE))

      if(any(grepl(":", f.notrans))) sd.x <- c(sd.x, sdInt(j, newdata.))

      sd.y <- sdFam(f.notrans[1], j, newdata.)

      if(intercepts == FALSE) Bnew <- B * (sd.x / sd.y) else

        Bnew <- c(0, B[-1] * (sd.x / sd.y))

    }

  } ) )

  ret <- data.frame(ret, Std.Estimate = unname(Bnew))

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
sdFam <- function(x, model, newdata) {

  .family = try(family(model), silent = TRUE)

  if(class(.family) == "try-error" & any(class(model) %in% c("glmerMod", "glmmPQL")))

     .family <-model$family

  if(class(.family) == "try-error" & !any(class(model) %in% c("glmerMod", "glmmPQL")))

    y <- newdata[, x] else

      if(.family$family == "gaussian") y <- newdata[, x] else

        y <- NA

  sd(y, na.rm = TRUE)

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

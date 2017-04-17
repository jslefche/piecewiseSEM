#' R^2 values from generalized linear mixed models
#'
#' @params a model or list of models

rsquared <- function(modelList, method = "trigamma") {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList = list(modelList)

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "formula", "formula.cerror")))]

  evaluateClasses(modelList)

  ret <- do.call(rbind, lapply(modelList, function(i) {

    if(all(class(i) %in% c("lm"))) r <- rsquared.lm(i)

    if(all(class(i) %in% c("gls"))) r <- rsquared.gls(i)

    if(any(class(i) %in% c("glm"))) r <- rsquared.glm(i)

    if(all(class(i) %in% c("lme"))) r <- rsquared.lme(i)

    if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) r <- rsquared.merMod(i)

    if(any(class(i) %in% c("glmerMod"))) r <- rsquared.glmerMod(i, method)

    if(any(class(i) %in% c("glmmPQL"))) r <- rsquared.glmmPQL(i, method)

    do.call(data.frame, r)

  } ) )

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 2)

  return(ret)

}

#' R^2 for lm objects
rsquared.lm <- function(model)

  list(family = "gaussian", link = "identity", Rsquared = summary(model)$r.squared)

#' R^2 for gls objects
rsquared.gls <- function(model) {

  X <- model.matrix(eval(model$call$model[-2]), getData.(model))

  sigmaF <- var(as.vector(model$coefficients %*% t(X)))

  sigmaE <- var(resid(model))

  list(family = "gaussian", link = "identity", Rsquared = sigmaF / (sigmaF + sigmaE))

}

#' R^2 for glm objects
rsquared.glm <- function(model) {

  if(family(model) == "gaussian") rsquared.lm(model) else {

    # Nagelkerke

    logLikModel <- as.numeric(logLik(model))

    null <- update(model, . ~ 1)

    logLikNull <- as.numeric(logLik(null))

    n <- 2/nobs(model)

    r <- 1 - (logLikNull/logLikModel)^n

    # why isnt this working!!

    r / (1 - logLikNull^n)

  }

}

#' R^2 for merMod objects
rsquared.merMod <- function(model) {

  X <- model.matrix(model)

  sigmaF <- var(as.vector(lme4::fixef(model) %*% t(X)))

  sigma <- unclass(lme4::VarCorr(model))

  sigmaL <- sum(sapply(1:length(sigma), function(i) {

    sigma. <- sigma[[i]]

    Z <- as.matrix(X[, rownames(sigma.), drop = FALSE])

    sum(diag(Z %*% sigma. %*% t(Z))) / nrow(X)

  } ) )

  sigmaE <- attr(lme4::VarCorr(model), "sc")^2

  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)

  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)

  list(family = "gaussian", link = "identity", Marginal = mar, Conditional = con)

}

#' R^2 for lme objects
rsquared.lme <- function(model) {

  X <- model.matrix(eval(model$call$fixed)[-2], model$data)

  sigmaF <- var(as.vector(nlme::fixef(model) %*% t(X)))

  # does not work for nested random effects! see maybe mgcv::extract.lme.cov ?

  sigma <- nlme::getVarCov(model)

  if(any(class(sigma) != "list")) sigma <- list(sigma)

  sigmaL <- sum(sapply(1:length(sigma), function(i) {

    sigma. <- sigma[[i]]

    Z <- as.matrix(X[, rownames(sigma.), drop = FALSE])

    sum(diag(Z %*% sigma. %*% t(Z))) / nrow(X)

  } ) )

  sigmaE <- summary(model)$sigma

  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)

  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)

  list(family = "gaussian", link = "identity", Marginal = mar, Conditional = con)

}

#' R^2 for glmer objects
rsquared.glmerMod <- function(model, method = "trigamma") {

  link <- family(model)$link

  family. <- family(model)$family

  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)

  if(family. == "Negative Binomial") rsquared.negbin(model, method) else {

    X <- model.matrix(model)

    sigmaF <- var(as.vector(lme4::fixef(model) %*% t(X)))

    if(family. == "poisson") {

      if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

      lambda <- attr(VarCorr(model), "sc")^2

      omega <- 1

      if(link == "mu^0.5") sigmaE <- 0.25 * omega else {

        if(link == "log") {

          nu <- omega / lambda

          if(method == "delta") sigmaE <- nu

          if(method == "lognormal") sigmaE <- log(1 + nu)

          if(method == "trigamma") sigmaE <- trigamma(nu)

        } else stop("Unsupported link function!")

      }

    } else

      if(family. == "binomial") {

        if(!method %in% c("none", "delta", "trigamma")) stop("Unsupported method!")

        lambda <- mean(model@resp$y)

        if(method == "trigamma") method = "delta"

        if(method == "none") sigmaE <- sigmaD <- (pi^2)/3

        if(method == "delta") sigmaE <- 1 / (lambda * (1 - lambda))

        } else

          if(family. == "Gamma") {

            if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

            lambda <- attr(VarCorr(model), "sc")^2

            omega <- 1

            if(link == "log") {

              nu <- omega / lambda

              if(method == "delta") sigmaE <- 1/nu

              if(method == "lognormal") sigmaE <- log(1 + 1/nu)

              if(method == "trigamma") sigmaE <- trigamma(nu)

              } else

                if(link == "inverse") {

                  } else stop("Unsupported link function!")

          } else stop("Unsupported family!")

    }

  sigmaA <- sum(as.numeric(VarCorr(model)))

  mar <- (sigmaF) / (sigmaF + sigmaA + sigmaE)

  con <- (sigmaF + sigmaA) / (sigmaF + sigmaA + sigmaE)

  list(family = family., link = link, Marginal = mar, Conditional = con)

}

#' R^2 for negbin objects
rsquared.negbin <- function(model, method = "trigamma") {

  link <- family(model)$link

  family. <- family(model)$family

  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)

  X <- model.matrix(model)

  sigmaF <- var(as.vector(fixef(model) %*% t(X)))

  if(family. == "Negative Binomial") {

    theta <- getME(model, "glmer.nb.theta")

    lambda <- mean(model@resp$y)

    if(link == "log") {

      if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

      nu <- (1/lambda) + (1/theta)

      if(method == "delta") sigmaE <- nu

      if(method == "lognormal") sigmaE <- log(1 + nu)

      if(method == "trigamma") sigmaE <- trigamma(nu^(-1))

    }

  }

  sigmaA <- sum(as.numeric(VarCorr(model)))

  mar <- (sigmaF) / (sigmaF + sigmaA + sigmaE)

  con <- (sigmaF + sigmaA) / (sigmaF + sigmaA + sigmaE)

  list(family = family., link = link, Marginal = mar, Conditional = con)

}

#' R^2 for glmmPQL objects
rsquared.glmmPQL <- function(model, method = "trigamma") {

  link <- model$family$link

  family. <- model$family$family

  X <- model.matrix(eval(model$call$fixed)[-2], model$data)

  sigmaF <- var(as.vector(nlme::fixef(model) %*% t(X)))

  sigma <- VarCorr(model)

  if(family. %in% c("poisson", "quasipoisson")) {

    lambda <- mean(model$data[, all.vars.merMod(formula(model))[1]])

    if(family. == "poisson") omega <- 1 else omega <- as.numeric(sigma[nrow(sigma), 1])

    if(link == "mu^0.5") sigmaE <- 0.25 * omega else {

      if(link == "log") {

        if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

        if(method == "delta") sigmaE <- omega / lambda

        if(method == "lognormal") sigmaE <- log(1 + (omega / lambda))

        if(method == "trigamma") sigmaE <- trigamma(lambda/omega)

      } else stop("Unsupported link function!")

    }

  } else if(family. %in% c("binomial", "quasibinomial")) {

    if(!method %in% c("none", "delta", "trigamma")) stop("Unsupported method!")

    lamba <- mean(model$data[, all.vars.merMod(formula(model))[1]])

    if(method == "trigamma") method = "delta"

    if(method == "none") sigmaE <- sigmaD <- pi^(2/3)

    if(method == "delta") sigmaE <- 1 / (lamba * (1 - lamba))

    } else if(family. %in% c("Gamma")) {

      if(link == "log") {

        if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

        nu <- 1 / as.numeric(VarCorr(model)[nrow(VarCorr(model)), 1])

        if(method == "delta") sigmaE <- 1 / nu

        if(method == "lognormal") sigmaE <- log(1 + 1/nu)

        if(method == "trigamma") sigmaE <- trigamma(nu)

        } else stop("Unsupported link function!")

      } else stop("Unsupported family!")

  sigmaA <- sum(as.numeric(sigma[!grepl("=|Residual", rownames(sigma)), 1]))

  mar <- (sigmaF) / (sigmaF + sigmaA + sigmaE)

  con <- (sigmaF + sigmaA) / (sigmaF + sigmaA + sigmaE)

  list(family = family., link = link, Marginal = mar, Conditional = con)

}

#' R^2 for glmmadmb objects




### DELETE LATER ###

#' Remove random effects from all.vars
all.vars.merMod <- function(.formula) {

  if(class(.formula) == "formula.cerror")

    gsub(" " , "", unlist(strsplit(.formula, "~~"))) else {

      n <- rownames(attr(terms(.formula), "factors"))

      if(any(grepl("\\|", n))) {

        f <- lme4::nobars(.formula)

        all.vars(f)

      } else all.vars(.formula)

    }

}

#' Get data from model list
getData. <- function(modelList) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  modelList <- modelList[!sapply(modelList, function(x) any(class(x) %in% c("matrix", "data.frame", "formula", "formula.cerror")))]

  data <- do.call(cbind, lapply(modelList, function(model) {

    if(any(class(model) %in% c("lm")))

      data <- model$model

    if(any(class(model) %in% c("glm", "glmmPQL")))

      data <- model$data

    if(any(class(model) %in% c("gls", "lme")))

      data <- nlme::getData(model)

    if(any(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod")))

      data <- model@frame

    return(data)

  } ) )

  data <- data[, !duplicated(colnames(data), fromLast = TRUE)]

  return(data)

}

#' Evaluate model classes and stop if unsupported model class
evaluateClasses <- function(modelList) {

  classes <- unlist(sapply(modelList, class))

  classes <- classes[!duplicated(classes)]

  model.classes <- c(
    "data.frame",
    "formula", "formula.cerror",
    "lm", "glm", "gls",
    "lme", "glmmPQL",
    "lmerMod", "merModLmerTest", "glmerMod"
  )

  if(!all(classes %in% model.classes))

    stop(
      paste0(
        "Unsupported model class in model list: ",
        paste0(classes[!classes %in% model.classes], collapse = ", "),
        ". See 'help(piecewiseSEM)' for more details."),
      call. = FALSE
    )

}

#' R-squared for linear regression
#'
#' Returns (pseudo)-R^2 values for all linear, generalized linear, and
#' generalized linear mixed effects models.
#'
#' For mixed models, marginal R2 considers only the variance by the fixed
#' effects, and the conditional R2 by both the fixed and random effects.
#'
#' For GLMs (\code{glm}), supported methods include: \itemize{
#' \item\code{mcfadden} 1 - ratio of likelihoods of full vs. null models
#'
#' \item\code{coxsnell} McFadden's R2 but raised to 2/N. Upper limit is < 1
#'
#' \item\code{nagelkerke} Adjusts Cox-Snell R2 so that upper limit = 1. The
#' DEFAULT method
#'
#' } For GLMERs fit to Poisson, Gamma, and negative binomial distributions
#' (\code{glmer}, \code{glmmPQL}, \code{glmer.nb}), supported methods include
#' \itemize{ \item\code{delta} Approximates the observation variance based on
#' second-order Taylor series expansion. Can be used with many families and
#' link functions
#'
#' \item\code{lognormal} Observation variance is the variance of the log-normal
#' distribution
#'
#' \item\code{trigamma} Provides most accurate estimate of the observation
#' variance but is limited to only the log link. The DEFAULT method
#'
#' } For GLMERs fit to the binomial distribution (\code{glmer},
#' \code{glmmPQL}), supported methods include: \itemize{
#' \item\code{theoretical} Assumes observation variance is pi^2/3
#'
#' \item\code{delta} Approximates the observation variance as above. The
#' DEFAULT method
#'
#' }
#'
#' @param modelList a regression, or a list of structural equations.
#' @param method The method used to compute the R2 value (See Details)
#' 
#' @return Returns a \code{data.frame} with the response, its family and link,
#' the method used to estimate R2, and the R2 value itself. Mixed models also
#' return marginal and conditional R2 values.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @references Nakagawa, Sinichi, Johnson, Paul C.D., and Holger Schielzeth.
#' "The coefficient of determination R2 and intra-class correlation coefficient
#' from generalized linear mixed-effects models revisted and expanded." bioRxiv
#' 095851 (2017).
#' @examples
#'
#'   \dontrun{
#'     # Create data
#'     dat <- data.frame(
#'       ynorm = rnorm(100),
#'       ypois = rpois(100, 100),
#'       x1 = rnorm(100),
#'       random = letters[1:5]
#'     )
#'
#'     # Get R2 for linear model
#'     rsquared(lm(ynorm ~ x1, dat))
#'
#'     # Get R2 for generalized linear model
#'     rsquared(glm(ypois ~ x1, "poisson", dat))
#'
#'     rsquared(glm(ypois ~ x1, "poisson", dat), method = "mcfadden") # McFadden R2
#'
#'     # Get R2 for generalized least-squares model
#'     rsquared(gls(ynorm ~ x1, dat))
#'
#'     # Get R2 for linear mixed effects model (nlme)
#'     rsquared(nlme::lme(ynorm ~ x1, random = ~ 1 | random, dat))
#'
#'     # Get R2 for linear mixed effects model (lme4)
#'     rsquared(lme4::lmer(ynorm ~ x1 + (1 | random), dat))
#'
#'     # Get R2 for generalized linear mixed effects model (lme4)
#'     rsquared(lme4::glmer(ypois ~ x1 + (1 | random), family = poisson, dat))
#'
#'     rsquared(lme4::glmer(ypois ~ x1 + (1 | random), family = poisson, dat), method = "delta")
#'
#'     # Get R2 for generalized linear mixed effects model (glmmPQL)
#'     rsquared(MASS::glmmPQL(ypois ~ x1, random = ~ 1 | random, family = poisson, dat))
#'   }
#'
#' @export
#' 
rsquared <- function(modelList, method = NULL) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList = list(modelList)

  modelList <- removeData(modelList, formulas = 1)

  evaluateClasses(modelList)

  ret <- lapply(modelList, function(i) {

    if(all(class(i) %in% c("lm", "pgls"))) r <- rsquared.lm(i) else

    if(all(class(i) %in% c("gls"))) r <- rsquared.gls(i) else

    if(any(class(i) %in% c("glm"))) r <- rsquared.glm(i, method) else

    # if(any(class(i) %in% c("phylolm", "phyloglm"))) r <- rsquared.phylolm(i)

    if(all(class(i) %in% c("lme"))) r <- rsquared.lme(i) else

    if(all(class(i) %in% c("lmerMod", "merModLmerTest", "lmerModLmerTest"))) r <- rsquared.merMod(i) else

    if(any(class(i) %in% c("glmerMod"))) r <- rsquared.glmerMod(i, method) else

    if(any(class(i) %in% c("glmmPQL"))) r <- rsquared.glmmPQL(i, method) else
      
      r <- list(family = "gaussian", link = "identity", method = "none", R.squared = NA)

    ret <- do.call(data.frame, r)

    ret <- data.frame(Response = all.vars.merMod(formula(i))[1], ret)

    # if(ncol(ret) != 5) ret[, ncol(ret) + 1] <- NA

    return(ret)

  } )

  if(length(unique(sapply(ret, ncol))) != 1) {

    nc <- max(sapply(ret, ncol))

    ret <- lapply(ret, function(i)

      if(ncol(i) == nc) i else {

        data.frame(
          Response = i$Response,
          family = i$family,
          link = i$link,
          method = i$method,
          Marginal = i$R.squared,
          Conditional = NA
        )

      }

    )

  }

  ret <- do.call(rbind, ret)

  # ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 2)

  rownames(ret) <- NULL

  return(ret)

}

#' R^2 for lm objects
#' 
#' @keywords internal
#' 
rsquared.lm <- function(model)

  list(family = "gaussian", link = "identity", method = "none", R.squared = summary(model)$r.squared)

#' R^2 for gls objects
#' 
#' @keywords internal
#' 
rsquared.gls <- function(model) {

  X <- model.matrix(eval(model$call$model[-2]), GetData(model))

  sigmaF <- var(as.vector(model$coefficients %*% t(X)))

  sigmaE <- var(resid(model))

  list(family = "gaussian", link = "identity", method = "none", R.squared = sigmaF / (sigmaF + sigmaE))

}

#' R^2 for glm objects
#' 
#' @keywords internal
#' 
rsquared.glm <- function(model, method = "nagelkerke") {

  if(is.null(method)) method <- "nagelkerke"

  link <- family(model)$link

  family. <- family(model)$family

  if(method %in% c("coxsnell", "nagelkerke")) {

    n <- nobs(model)

    r <- 1 - exp((model$dev - model$null) / n)

    if(method == "nagelkerke")

      r <- r / (1 - exp(-model$null / n))

    }

  if(method == "mcfadden") {

    null <- update(model, . ~ 1)

    r <- 1 - (as.numeric(logLik(model)) / as.numeric(logLik(null)))

  }

  list(family = family., link = link, method = method, R.squared = r)

}

#' R^2 for phylolm objects
#' 
#' @keywords internal
#' 
# rsquared.phylolm <- function(model) {
# 
#   family. <- ifelse(class(model) == "phylolm", "gaussian", model$method)
# 
#   link <- ifelse(class(model) == "phylolm", "identity", "?")
# 
#   list(family = family., link = NA, method = "none", R.squared = NA)
# 
# }

#' R^2 for merMod objects
#' 
#' @keywords internal
#' 
rsquared.merMod <- function(model) {

  X <- model.matrix(model)

  sigmaF <- var(as.vector(lme4::fixef(model) %*% t(X)))

  sigma <- unclass(lme4::VarCorr(model))

  sigmaL <- sum(sapply(1:length(sigma), function(i) {

    sigma. <- sigma[[i]]

    Z <- as.matrix(X[, rownames(sigma.), drop = FALSE])

    sum(rowSums(Z %*% sigma.) * Z) / nrow(X)

  } ) )

  sigmaE <- attr(lme4::VarCorr(model), "sc")^2

  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)

  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)

  list(family = "gaussian", link = "identity", method = "none", Marginal = mar, Conditional = con)

}

#' R^2 for lme objects
#' 
#' @keywords internal
#' 
rsquared.lme <- function(model) {

  X <- model.matrix(eval(model$call$fixed)[-2], droplevels(model$data))

  sigmaF <- var(as.vector(nlme::fixef(model) %*% t(X)))

  sigma <- GetVarCov(model)

  sigmaL <- sum(sapply(sigma, function(i) {

    Z <- as.matrix(X[, rownames(i), drop = FALSE])

    sum(rowSums(Z %*% i) * Z) / nrow(X)

  } ) )

  sigmaE <- summary(model)$sigma^2

  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)

  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)

  list(family = "gaussian", link = "identity", method = "none", Marginal = mar, Conditional = con)

}

#' R^2 for glmer objects
#' 
#' @keywords internal
#' 
rsquared.glmerMod <- function(model, method = "trigamma") {

  if(is.null(method)) method <- "trigamma"

  link <- family(model)$link

  family. <- family(model)$family

  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)

  if(family. == "Negative Binomial") rsquared.negbin(model, method) else {

    X <- model.matrix(model)

    sigmaF <- var(as.vector(lme4::fixef(model) %*% t(X)))

    sigma <- unclass(lme4::VarCorr(model))

    data <- GetData(model)
    
    if(family. == "poisson") {

      if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

      sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "RE"))
      
      sigmaE <- sum(GetOLRE(sigma, model, X, data, RE = "OLRE"))
      
      rand <- onlyBars(formula(model))
      
      f <- paste(all.vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
      
      nullmodel <- suppressWarnings(lme4::glmer(formula(f), family = poisson(link = link), data = data))
      
      # lambda <- attr(VarCorr(model), "sc")^2
      
      lambda <- exp(fixef(nullmodel)[1] + (sigmaL + sigmaE)/2)

      omega <- 1

      if(link == "mu^0.5") sigmaD <- 0.25 * omega else {

        if(link == "log") {

          nu <- omega / lambda

          if(method == "delta") sigmaD <- nu

          if(method == "lognormal") sigmaD <- log(1 + nu)

          if(method == "trigamma") sigmaD <- trigamma(1/nu)

        } else stop("Unsupported link function!")

      }

    } else

      if(family. == "binomial") {

        if(method == "trigamma") method <- "delta"

        if(!method %in% c("theoretical", "delta")) stop("Unsupported method!")
        
        sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "all"))
        
        sigmaE <- 0

        if(method == "theoretical") sigmaD <- pi^2/3

        if(method == "delta") {
          
          rand <- onlyBars(formula(model))
          
          f <- paste(all.vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
          
          nullmodel <- suppressWarnings(lme4::glmer(formula(f), family = binomial(link = link), data = data))
          
          vt <- sum(unlist(VarCorr(nullmodel)))

          pmean <- plogis(as.numeric(fixef(nullmodel)) - 0.5 * vt *
                            tanh(as.numeric(fixef(nullmodel)) * (1 + 2 * exp(-0.5 * vt))/6))

          sigmaD <- 1/(pmean * (1 - pmean))

          }

        } else

          if(family. == "Gamma") {

            if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

            sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "all"))
            
            sigmaE <- 0
            
            lambda <- attr(VarCorr(model), "sc")^2

            omega <- 1

            if(link == "log") {

              nu <- omega / lambda

              if(method == "delta") sigmaD <- 1/nu

              if(method == "lognormal") sigmaD <- log(1 + 1/nu)

              if(method == "trigamma") sigmaD <- trigamma(nu)

              } else

                if(link == "inverse") {
                  
                  if(method != "delta") method <- "delta"
                  
                  nu <- omega / lambda
                  
                  sigmaD <- 1/(nu * lambda^2)

                  } else stop("Unsupported link function!")

          } else stop("Unsupported family!")

    mar <- (sigmaF) / (sigmaF + sigmaL + sigmaD + sigmaE)
    
    con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaD + sigmaE)

    list(family = family., link = link, method = method, Marginal = mar, Conditional = con)

    }

}

#' R^2 for negbin objects
#' 
#' @keywords internal
#' 
rsquared.negbin <- function(model, method = "trigamma") {

  if(is.null(method)) method <- "trigamma"

  link <- family(model)$link

  family. <- family(model)$family

  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)

  X <- model.matrix(model)

  sigmaF <- var(as.vector(fixef(model) %*% t(X)))

  sigma <- unclass(lme4::VarCorr(model))

  rand <- sapply(lme4::findbars(formula(model)), function(x) as.character(x)[3])

  rand <- rand[!duplicated(rand)]

  data <- GetData(model)

  idx <- sapply(strsplit(rand, "\\:"), function(x) {

    length(unique(data[, x])) == nrow(data)

  } )

  sigma.names <- strsplit(names(sigma), "\\.")

  idx. <- sapply(sigma.names, function(x) !any(x %in% rand[idx]))

  sigmaL <- sum(sapply(sigma[idx.], function(i) {

    Z <- as.matrix(X[, rownames(i), drop = FALSE])

    sum(rowSums(Z %*% i) * Z) / nrow(X)

  } ) )

  if(family. == "Negative Binomial") {

    rand <- onlyBars(formula(model))

    f <- paste(all.vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))

    nullmodel <- suppressWarnings(lme4::glmer.nb(formula(f), family=negative.binomial(link = link), data))

    # nullmodel <- update(model, formula(paste(". ~ 1 +", onlyBars(formula(model)))))

    lambda <- as.numeric(exp(fixef(nullmodel) + 0.5 * sum(unlist(VarCorr(nullmodel)))))

    theta <- lme4::getME(model, "glmer.nb.theta")

    if(link == "log") {

      if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

      nu <- (1/lambda) + (1/theta)

      if(method == "delta") sigmaE <- nu

      if(method == "lognormal") sigmaE <- log(1 + nu)

      if(method == "trigamma") sigmaE <- trigamma(nu^(-1))

    }

  }

  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)

  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)

  list(family = family., link = link, method = method, Marginal = mar, Conditional = con)

}

#' R^2 for glmmPQL objects
#' 
#' @keywords internal
#' 
rsquared.glmmPQL <- function(model, method = "trigamma") {

  if(is.null(method)) method <- "trigamma"

  link <- model$family$link

  family. <- model$family$family

  X <- model.matrix(eval(model$call$fixed)[-2], droplevels(model$data))

  sigmaF <- var(as.vector(nlme::fixef(model) %*% t(X)))

  sigma <- VarCorr(model)

  sigmaL <- sum(as.numeric(sigma[!grepl("=|Residual", rownames(sigma)), 1]))

  data <- GetData(model)

  if(family. %in% c("poisson", "quasipoisson")) {

    f <- paste(all.vars_trans(formula(model))[1], " ~ 1")

    if(family. == "poisson")

      nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = poisson(link = link), data = data, verbose = FALSE)) else

        nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = quasipoisson(link = link), data = data, verbose = FALSE))

    lambda <- as.numeric(exp(fixef(nullmodel) + 0.5 * sum(unlist(GetVarCov(nullmodel)))))

    if(family. == "poisson") omega <- 1 else omega <- as.numeric(sigma[nrow(sigma), 1])

    if(link == "mu^0.5") sigmaE <- 0.25 * omega else {

      if(link == "log") {

        if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

        nu <- omega / lambda
        
        if(method == "delta") sigmaE <- nu

        if(method == "lognormal") sigmaE <- log(1 + (nu))

        if(method == "trigamma") sigmaE <- trigamma(1/nu)

      } else stop("Unsupported link function!")

    }

  } else if(family. %in% c("binomial", "quasibinomial")) {

    if(method == "trigamma") method <- "delta"

    if(!method %in% c("theoretical", "delta")) stop("Unsupported method!")

    if(method == "theoretical") sigmaE <- sigmaD <- pi^2/3

    if(method == "delta") {

      f <- paste(all.vars_trans(formula(model))[1], " ~ 1")

      if(family. == "binomial")

        nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = binomial(link = link), data = data, verbose = FALSE)) else

          nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = quasibinomial(link = link), data = data, verbose = FALSE))

      vt <- sum(unlist(GetVarCov(nullmodel)))

      pmean <- plogis(as.numeric(fixef(nullmodel)) - 0.5 * vt *
                        tanh(as.numeric(fixef(nullmodel)) * (1 + 2 * exp(-0.5 * vt))/6))

      sigmaE <- 1/(pmean * (1 - pmean))

      }

    } else if(family. %in% c("Gamma")) {

      if(link == "log") {

        if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")

        nu <- 1 / as.numeric(VarCorr(model)[nrow(VarCorr(model)), 1])

        if(method == "delta") sigmaE <- 1 / nu

        if(method == "lognormal") sigmaE <- log(1 + 1/nu)

        if(method == "trigamma") sigmaE <- trigamma(nu)

        } else stop("Unsupported link function!")

      } else stop("Unsupported family!")

  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)

  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)

  list(family = family., link = link, method = method, Marginal = mar, Conditional = con)

}

# R^2 for glmmadmb objects
# 
# @keywords internal
#


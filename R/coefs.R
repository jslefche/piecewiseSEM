#' Extract path coefficients
#'
#' Extracts (standardized) path coefficients from a \code{psem} object.
#'
#' P-values for models constructed using \code{lme4} are obtained
#' using the Kenward-Roger approximation of the denominator degrees of freedom
#' as implemented in the \code{pbkrtest} package.
#'
#' Different forms of standardization can be implemented using the \code{standardize}
#' argument:\itemize{
#' \item{\code{none} No standardized coefficients are reported.}
#' \item{\code{scale} Raw coefficients are scaled by the ratio of the standard deviation
#' of x divided by the standard deviation of y. See below for cases pertaining to GLM. }
#' \item{\code{range} Raw coefficients are scaled by a pre-selected range of x
#' divided by a preselected range of y. The default argument is \code{range} which takes the
#' two extremes of the data, otherwise the user must supply must a named \code{list} where
#' the names are the variables to be standardized, and each entry contains a vector of length 2
#' to the ranges to be used in standardization.}
#' }
#'
#' For binary response models (i.e., binomial responses), standardized coefficients
#' are obtained in one of two ways:\itemize{
#' \item{\code{latent.linear} Referred to in Grace et al. (in review) as the standard form of
#' the latent-theoretic (LT) approach. In this method, there is assumed to be a continuous
#' latent propensity, y*, that underlies the observed binary responses. The standard
#' deviation of y* is computed as the square-root of the variance of the predictions
#' (on the linear or 'link' scale) plus the distribution-specific assumed variance
#' (for logit links: pi^2/3, for probit links: 1).}
#' \item{\code{Menard.OE} Referred to in Grace et al. (in review) as the standard form of
#' the observed-empirical (OE) approach. In this method, error variance is based on the
#' differences between predicted scores and the observed binary data. The standard
#' deviation used for standardization is computed as the square-root of the variance of
#' the predictions (on the linear scale) plus the correlation between the observed and
#' predicted (on the original or 'response' scale) values of y.}
#' }
#'
#' @param modelList A list of structural equations.
#' @param standardize The type of standardization: \code{none}, \code{scale}, \code{range}.
#' Default is \code{scale}.
#' @param standardize.type The type of standardized for non-Gaussian responses:
#' \code{latent.linear}, \code{Menard.OE}. Default is \code{latent.linear}.
#' @param intercepts Whether intercepts should be included in the coefficients
#' table. Default is FALSE.
#' @return Returns a \code{data.frame} of coefficients, their standard errors,
#' degrees of freedom, and significance tests.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>, Jim Grace
#' @references Grace, J.B., Johnson, D.A., Lefcheck, J.S., and Byrnes, J.E.
#' "Standardized Coefficients in Regression and Structural Models with Binary Outcomes." 
#' Ecosphere 9(6): e02283.
#' @seealso \code{\link{KRmodcomp}}
#' 
#' @export
#'
coefs <- function(modelList, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(class(modelList) == "psem") data <- modelList$data else data <- GetData(modelList)

  if(class(data) %in% c("SpatialPointsDataFrame")) data <- data@data

  if(class(data) %in% c("comparative.data")) data <- data$data

  if(all(standardize != "none")) ret <- stdCoefs(modelList, data, standardize, standardize.type, intercepts) else

    ret <- unstdCoefs(modelList, data, intercepts)

  ret <- cbind(ret, isSig(ret$P.Value))

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)

  names(ret)[length(ret)] <- ""

  return(ret)

}

#' Get raw (understandardized) coefficients from model
#' 
#' @keywords internal
#' 
#' @export
#' 
unstdCoefs <- function(modelList, data = NULL, intercepts = FALSE) {
  
  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)
  
  if(is.null(data) & class(modelList) == "psem") data <- modelList$data
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 2)
  
  ret <- do.call(rbind, lapply(modelList, function(i) {
    
    if(all(class(i) %in% c("formula.cerror")))
      
      ret <- cerror(i, modelList, data) else {
        
        if(all(class(i) %in% c("lm", "glm", "negbin", "lmerMod", "glmerMod", "merModLmerTest", "pgls", "phylolm", "phyloglm"))) {
          
          ret <- as.data.frame(summary(i)$coefficients)
          
          if(all(class(i) %in% c("lm", "glm", "negbin"))) ret <- cbind(ret[, 1:2], DF = summary(i)$df[2], ret[, 3:4])
          
          if(all(class(i) %in% c("glmerMod", "pgls"))) ret <- cbind(ret[, 1:2], DF = length(summary(i)$residuals), ret[, 3:4])
          
          if(all(class(i) %in% c("phylolm", "phyloglm"))) ret <- cbind(ret[, 1:2], DF = i$n, ret[, c(3, 6)])
          
          if(all(class(i) %in% c("lmerMod", "merModLmerTest"))) {
            
            krp <- KRp(i, all.vars_trans(formula(i))[-1], data, intercepts = intercepts)
            
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
          Response = all.vars_trans(listFormula(list(i))[[1]])[1],
          Predictor = rownames(ret),
          ret
        )
        
        names(ret) <- c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value")
        
      }
    
    if(intercepts == FALSE) ret <- ret[ret$Predictor != "(Intercept)", ]
    
    return(ret)
    
  } ) )
  
  rownames(ret) <- NULL
  
  return(ret)
  
}

#' Calculate standardized regression coefficients
#' 
#' @keywords internal
#' 
#' @export
#' 
stdCoefs <- function(modelList, data = NULL, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE) {
  
  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)
  
  if(is.null(data) & class(modelList) == "psem") data <- modelList$data
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 2)
  
  ret <- unstdCoefs(modelList, data, intercepts)
  
  Bnew <- do.call(rbind, lapply(1:length(modelList), function(i) {
    
    j <- modelList[[i]]
    
    f <- listFormula(list(j))[[1]]
    
    newdata <- data[, all.vars.merMod(f)]
    
    f.trans <- all.vars_trans(f)
    
    f.notrans <- all.vars_notrans(f)
    
    if(all(class(j) %in% c("formula.cerror"))) {
      
      data.frame(Std.Estimate = ret[ret$Response == paste0("~~", f.trans[1]) & 
                                      ret$Predictor == paste0("~~", f.trans[2]), "Estimate"])
      
    } else {
      
      if(any(class(newdata) %in% c("SpatialPointsDataFrame"))) newdata <- newdata@data
      
      newdata <- dataTrans(formula(j), newdata)
      
      B <- ret[ret$Response == f.trans[1], "Estimate"]
      
      # B.se <- ret[ret$Response == f.trans[1]), "Std.Error"]
      
      if(all(standardize == "scale"))
        
        sd.x <- sapply(f.notrans[!grepl(":", f.notrans)][-1], function(x) sd(newdata[, x], na.rm = TRUE)) else
          
          if(all(standardize == "range"))
            
            sd.x <- sapply(f.notrans[!grepl(":", f.notrans)][-1], function(x) diff(range(newdata[, x], na.rm = TRUE))) else
              
              if(is.list(standardize)) {
                
                vars <- unlist(sapply(modelList, all.vars_notrans))
                
                vars <- vars[!grepl(":", vars)]
                
                if(!all(names(standardize) %in% vars))
                  
                  stop("Names in standardize list must match those in the model formula!")
                
                sd.x <- sapply(f.notrans[!grepl(":", f.notrans)][-1], function(x) {
                  
                  nm <- which(names(standardize) == x)
                  
                  if(sum(nm) == 0) {
                    
                    warning(paste0("Relevant range not specified for variable '", x, "'. Using observed range instead"), call. = FALSE)
                    
                    diff(range(newdata[, x], na.rm = TRUE)) 
                    
                  } else  diff(range(standardize[[nm]]))
                  
                } )
                
              } else stop("`standardize` must be either 'scale' or 'range' (or a list of ranges).", call. = FALSE)
      
      if(any(grepl(":", f.notrans))) sd.x <- c(sd.x, scaleInt(j, newdata, standardize))
      
      sd.y <- scaleFam(f.notrans[1], j, newdata, standardize, standardize.type)
      
      if(intercepts == FALSE)
        
        data.frame(Std.Estimate = B * (sd.x / sd.y)) else
          
          data.frame(Std.Estimate = c(0, B[-1] * (sd.x / sd.y)))
      
    }
    
  } ) )

ret <- data.frame(ret, Bnew)

rownames(ret) <- NULL

return(ret)

}

#' Transform variables based on model formula and store in new data frame
#' 
#' @keywords internal
#' 
dataTrans <- function(formula., newdata) {

  notrans <- all.vars.merMod(formula.)

  trans <- all.vars_trans(formula.)

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
#' 
#' @keywords internal
#' 
scaleFam <- function(y, model, newdata, standardize = "scale", standardize.type = "latent.linear") {

  family. <- try(family(model), silent = TRUE)

  if(class(family.) == "try-error") family. <- try(model$family, silent = TRUE)

  if(class(family.) == "try-error" | is.null(family.) & all(class(model) == "lme"))

    family. <- list(family = "gaussian", link = "identity")

  if(class(family.) == "try-error" | is.null(family.) | any(class(model) %in% c("glmerMod", "glmmPQL")))

    sd.y <- NA else {

      if(family.$family == "gaussian") {

        if(all(standardize == "scale")) sd.y <- sd(newdata[, y], na.rm = TRUE) else

          if(all(standardize == "range")) sd.y <- diff(range(newdata[, y], na.rm = TRUE)) else

            if(is.list(standardize)) {

              nm <- which(names(standardize) == y)
              
              if(sum(nm) == 0) {
                
                warning(paste0("Relevant range not specified for variable '", y, "'. Using observed range instead"), call. = FALSE)
                
                sd.y <- diff(range(newdata[, y], na.rm = TRUE)) 
                
              } else sd.y <- diff(range(standardize[[nm]]))
              
            }

          } else if(family.$family == "binomial")

              sd.y <- scaleGLM(model, standardize, standardize.type) else

                sd.y <- NA

    }

  return(sd.y)

}

#' Compute standard deviation or relevant range of response for GLMs
#' 
#' @keywords internal
#' 
scaleGLM <- function(model, standardize = "scale", standardize.type = "latent.linear") {

  preds <- predict(model, type = "link")

  if(standardize.type == "Menard.OE") {

    y <- all.vars_notrans(model)[1]

    data <- getSingleData(model)

    R <- cor(data[, y], predict(model, type = "response"))

    sd.y <- sqrt(var(preds)) / R

  }

  if(standardize.type == "latent.linear") {

    link. <- family(model)$link

    if(link. == "logit") sigmaE <- pi^2/3 else

      if(link. == "probit") sigmaE <- 1

    sd.y <- sqrt(var(preds) + sigmaE)

  }

  if(all(standardize == "range") | is.list(standardize)) sd.y <- sd.y * 6

  return(sd.y)

}

#' Calculate standard deviation or relevant range for interaction terms
#' 
#' @keywords internal
#' 
scaleInt <- function(model, newdata, standardize) {

  v <- attr(terms(model), "term.labels")

  int <- v[grepl(":", v)]

  sapply(int, function(x) {

    x <- strsplit(x, ":")[[1]]

    x <- gsub("(.*) \\+.*", "\\1", gsub(".*\\((.*)\\)", "\\1", x))

    p <- apply(newdata[, x], 1, prod, na.rm = TRUE)

    if(standardize == "scale") sd(p) else if(standardize == "range")

      diff(range(p, na.rm = TRUE)) else if(is.list(standardize))
        
        stop("Relevant range standardization not applicable to models with interactions!")

  } )

}

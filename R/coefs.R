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
#' the names are the variables to be standardized, and each entry contains a vector of 
#' length == 2 to the ranges to be used in standardization.}
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
#' For categorical predictors: significance is determined using ANOVA (or analysis of
#' deviance). Because n-1 coefficients are reported for n levels, the output instead
#' reports model-estimated means in the \code{Estimate} column. This is done so all
#' n paths in the corresponding path diagram have assignable values.
#' 
#' The means are generated using function \code{emmeans} in the \code{emmeans} package. 
#' Pairwise contrasts  are further conducted among all levels using the default  
#' correction for multiple testing. The results of those comparisons are given in the 
#' significance codes (e.g., "a", "b", "ab") as reported in the \code{emmeans::cld} function.
#'
#' @param modelList A list of structural equations, or a model.
#' @param standardize The type of standardization: \code{none}, \code{scale}, \code{range}.
#' Default is \code{scale}.
#' @param standardize.type The type of standardized for non-Gaussian responses:
#' \code{latent.linear}, \code{Menard.OE}. Default is \code{latent.linear}.
#' @param test.type the type of test for significance of categorical variables
#' from \code{\link{car::Anova}}. Default is type "III".
#' @param intercepts Whether intercepts should be included in the coefficients
#' table. Default is FALSE.
#' @return Returns a \code{data.frame} of coefficients, their standard errors,
#' degrees of freedom, and significance tests.
#' @author Jon Lefcheck <LefcheckJ@@si.edu>, Jim Grace
#' @references Grace, J.B., Johnson, D.A., Lefcheck, J.S., and Byrnes, J.E.
#' "Standardized Coefficients in Regression and Structural Models with Binary Outcomes." 
#' Ecosphere 9(6): e02283.
#' @seealso \code{\link{KRmodcomp}}, \code{\link{emmeans}}, \code{\link{CLD}}
#' @examples 
#' mod <- psem(
#' lm(rich ~ cover, data = keeley),
#' lm(cover ~ firesev, data = keeley),
#' lm(firesev ~ age, data = keeley),
#' data = keeley
#' )
#' 
#' coefs(mod)
#' 
#' @export
#'
coefs <- function(modelList, standardize = "scale", standardize.type = "latent.linear", test.type = "III", intercepts = FALSE) {

  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)

  if(class(modelList) == "psem") data <- modelList$data else data <- GetData(modelList)

  if(class(data) %in% c("SpatialPointsDataFrame")) data <- data@data

  if(class(data) %in% c("comparative.data")) data <- data$data

  if(all(standardize != "none")) ret <- stdCoefs(modelList, data, standardize, standardize.type, test.type, intercepts) else

    ret <- unstdCoefs(modelList, data, test.type, intercepts)

  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)

  ret[is.na(ret)] <- "-"
  
  names(ret)[length(ret)] <- ""

  return(ret)

}

#' Get raw (undstandardized) coefficients from model
#' 
#' @keywords internal
#' 
#' @export
#' 
unstdCoefs <- function(modelList, data = NULL, test.type = "II", intercepts = FALSE) {
  
  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)
  
  if(is.null(data) & class(modelList) == "psem") data <- modelList$data
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 2)
  
  ret <- do.call(rbind, lapply(modelList, function(i) {
    
    if(all(class(i) %in% c("formula.cerror"))) {
      
      ret <- cbind.data.frame(cerror(i, modelList, data)) 
      
      names(ret) <- c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value", "")

      } else {
        
        ret <- getCoefficients(i, data, test.type)
        
        if(intercepts == FALSE) ret <- ret[ret$Predictor != "(Intercept)", ]
        
      }
    
    return(ret)
    
  } ) )
  
  rownames(ret) <- NULL
  
  return(ret)
  
}

#' Get coefficients from linear regression
#' 
#' @keywords internal
#' 
#' @export
#'
getCoefficients <- function(model, data, test.type = "II") {
  
  vars <- all.vars.merMod(model)
  
  factorVars <- vars[which(sapply(data[, vars], class) == "factor")]

  if(all(class(model) %in% c("lm", "glm", "negbin", "lmerMod", "glmerMod", "lmerModLmerTest", "pgls", "phylolm", "phyloglm"))) {
    
    ret <- as.data.frame(summary(model)$coefficients)[rownames(summary(model)$coefficients) %in% c("(Intercept)", vars), ]
    
    if(all(class(model) %in% c("lm", "glm", "negbin"))) ret <- cbind(ret[, 1:2], DF = summary(model)$df[2], ret[, 3:4])
    
    if(all(class(model) %in% c("glmerMod", "pgls"))) ret <- cbind(ret[, 1:2], DF = length(summary(model)$residuals), ret[, 3:4])
    
    if(all(class(model) %in% c("phylolm", "phyloglm"))) ret <- cbind(ret[, 1:2], DF = model$n, ret[, c(3, 6)])
  
    if(all(class(model) %in% c("lmerMod"))) {

      krp <- KRp(model, vars[-1], data, intercepts = TRUE)
      
      ret <- cbind.data.frame(ret[, 1:2], DF = krp[, 1], ret[, 3, drop = FALSE], krp[, 2])
        
      names(ret)[ncol(ret)] <- "Pr(>|t|)"
      
    }
    
  }
  
  if(all(class(model) %in% c("sarlm"))) {
    
    ret <- as.data.frame(summary(model)$Coef)
    
    ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])
    
  }
  
  if(all(class(model) %in% c("gls", "lme", "glmmPQL"))) {
    
    ret <- as.data.frame(summary(model)$tTable)[rownames(summary(model)$tTable) %in% c("(Intercept)", vars), ]
    
    if(ncol(ret) == 4 & all(class(model) %in% c("gls")))
      
      ret <- cbind(ret[, 1:2], DF = length(residuals(model)), ret[, 3:4])
    
  }
  
  ret <- cbind(ret, isSig(ret[, 5]))
  
  if(length(factorVars) > 0) {
    
    anovaTable <- as.data.frame(car::Anova(model, type = test.type))
  
    levs <- lapply(factorVars, function(j) suppressMessages(emmeans::emmeans(model, list(formula(paste(" ~", j))))))
  
    ret <- do.call(rbind, lapply(levs, function(j) {
      
       out <- as.data.frame(emmeans::CLD(j[[1]],  Letters = letters[1:26]))
     
       rownames(out) <- paste0(names(out)[1], "[", out[, 1], "] mean=")
       
       atab <- anovaTable[which(grepl(names(out)[1], rownames(anovaTable))), ]
       
       atab. <- data.frame(emmean = NA, SE = NA, DF = atab$Df, crit.value = atab[, 1], P.value = atab[, ncol(atab)], sig = isSig(atab[, ncol(atab)]))
       
       rownames(atab.) <- rownames(atab) 
       
       out <- rbind(
         atab.,
         data.frame(out[, 2:3], DF = out$df, crit.value = NA, P.value = NA, sig = out[, 7])
       )
       
       names(out) <- names(ret)
       
       ret <- rbind(ret, out)
       
       } ) )
    
  }
  
  ret <- data.frame(
    Response = all.vars_trans(listFormula(list(model))[[1]])[1],
    Predictor = rownames(ret),
    ret
  )
  
  names(ret) <- c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value", "")
  
  rownames(ret) <- NULL
  
  ret$Response <- as.character(ret$Response)
  
  ret[is.na(ret$P.Value), "Response"] <- ""
  
  return(ret)  
     
}

#' Calculate standardized regression coefficients
#' 
#' @keywords internal
#' 
#' @export
#' 
stdCoefs <- function(modelList, data = NULL, standardize = "scale", standardize.type = "latent.linear", test.type = "II", intercepts = FALSE) {
  
  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)
  
  if(is.null(data) & class(modelList) == "psem") data <- modelList$data
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 2)
  
  ret <- do.call(rbind, lapply(modelList, function(i) {
    
    if(all(class(i) %in% c("formula.cerror"))) {
      
      ret <- cerror(i, modelList)
      
      cbind.data.frame(ret[, 1:7], Std.Estimate = ret[, 3], sig = ret[, 8])
      
    } else {
      
      ret <- unstdCoefs(i, data, test.type, intercepts)
      
      vars <- all.vars.merMod(i)
      
      newdata <- data[, vars]
      
      if(any(class(newdata) %in% c("SpatialPointsDataFrame"))) newdata <- newdata@data
      
      newdata <- dataTrans(formula(i), newdata)
      
      numVars <- vars[which(sapply(data[, vars], class) != "factor")]
      
      B <- ret[ret$Predictor %in% numVars[-1], "Estimate"]
      
      names(B) <- ret[ret$Predictor %in% numVars[-1], "Predictor"]
      
      sd.x <- GetSDx(i, modelList, data, standardize)
      
      sd.y <- GetSDy(i, newdata, standardize, standardize.type)
      
      if(intercepts == FALSE)
        
        Std.Estimate <- B * (sd.x / sd.y) else
          
          Std.Estimate <- c(0, B[-1] * (sd.x / sd.y))
      
      ret[which(ret$Predictor %in% numVars[-1]), "Std.Estimate"] <- Std.Estimate
      
      ret <- cbind(ret[, 1:7], ret[, 9, drop = FALSE], sig = ret[, 8])
      
      return(ret)
      
    }
    
  } ) )  
  
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

#' Get standard deviation of predictor variables
#' 
#' @keywords internal
#' 
GetSDx <- function(model, modelList, data, standardize = "scale") {
  
  vars <- all.vars.merMod(model)
  
  numVars <- vars[which(sapply(data[, vars], class) != "factor")]
  
  if(all(standardize == "scale"))
    
    sd.x <- sapply(numVars[!grepl(":", numVars)][-1], function(x) sd(data[, x], na.rm = TRUE)) else
      
      if(all(standardize == "range"))
        
        sd.x <- sapply(numVars[!grepl(":", numVars)][-1], function(x) diff(range(data[, x], na.rm = TRUE))) else
          
          if(is.list(standardize)) {
            
            vars <- unlist(sapply(modelList, all.vars_notrans))
            
            vars <- vars[!grepl(":", vars)]
            
            if(!all(names(standardize) %in% vars))
              
              stop("Names in standardize list must match those in the model formula!")
            
            sd.x <- sapply(numVars[!grepl(":", numVars)][-1], function(x) {
              
              nm <- which(names(standardize) == x)
              
              if(sum(nm) == 0) {
                
                warning(paste0("Relevant range not specified for variable '", x, "'. Using observed range instead"), call. = FALSE)
                
                diff(range(data[, x], na.rm = TRUE)) 
                
              } else  diff(range(standardize[[nm]]))
              
            } )
            
          } else stop("`standardize` must be either 'scale' or 'range' (or a list of ranges).", call. = FALSE)
  
  if(any(grepl(":", numVars))) sd.x <- c(sd.x, scaleInt(i, data, standardize))
  
  if(length(sd.x) == 0) sd.x <- NA
  
  return(sd.x)
  
}

#' Properly scale standard deviations depending on the error distribution
#' 
#' @keywords internal
#' 
GetSDy <- function(model, data, standardize = "scale", standardize.type = "latent.linear") {
  
  vars <- all.vars.merMod(model)
  
  y <- vars[1]
  
  family. <- try(family(model), silent = TRUE)
  
  if(class(family.) == "try-error") family. <- try(model$family, silent = TRUE)
  
  if(class(family.) == "try-error" | is.null(family.) & all(class(model) %in% c("sarlm", "lme")))
    
    family. <- list(family = "gaussian", link = "identity")
  
  if(class(family.) == "try-error" | is.null(family.) | any(class(model) %in% c("glmerMod", "glmmPQL")))
    
    sd.y <- NA else {
      
      if(family.$family == "gaussian") {
        
        if(all(standardize == "scale")) sd.y <- sd(data[, y], na.rm = TRUE) else
          
          if(all(standardize == "range")) sd.y <- diff(range(data[, y], na.rm = TRUE)) else
            
            if(is.list(standardize)) {
              
              nm <- which(names(standardize) == y)
              
              if(sum(nm) == 0) {
                
                warning(paste0("Relevant range not specified for variable '", y, "'. Using observed range instead"), call. = FALSE)
                
                sd.y <- diff(range(data[, y], na.rm = TRUE)) 
                
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

    data <- GetSingleData(model)

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

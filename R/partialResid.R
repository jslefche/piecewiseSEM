#' Computing partial effects
#'
#' Extracts partial residuals from a model or \code{psem} object for a given
#' \code{x} and \code{y}.
#'
#' This function computes the partial residuals of \code{y ~ x + Z} in a
#' two-step procedure to remove the variation explained by \code{Z}: (1) remove
#' \code{x} from the equation and model \code{y ~ Z}, and (2) replace \code{y}
#' with \code{x} and model \code{x ~ Z}.
#'
#' @param formula.  A formula where the \code{lhs} is the response and the
#' \code{rhs} is the predictor whose partial effect is desired.
#' @param modelList A list of structural equations.
#' @param data A \code{data.frame} used to fit the equations.
#' @return Returns a \code{data.frame} of residuals of \code{y ~ Z} called
#' \code{yresids}, of \code{x ~ Z} called \code{xresids}.
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso \code{\link{cerror}}
#' @examples
#'
#' # Generate data
#' dat <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
#'
#' # Build model
#' model <- lm(y ~ x1 + x2, dat)
#'
#' # Compute partial residuals of y ~ x1
#' yresid <- resid(lm(y ~ x2, dat))
#'
#' xresid <- resid(lm(x1 ~ x2, dat))
#'
#' plot(yresid, xresid)
#'
#' # Use partialResid
#' presid <- partialResid(y ~ x1, model)
#'
#' plot(presid) # identical plot!
#'
#' @export partialResid
#' 
partialResid <- function(formula., modelList, data = NULL) {
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  data <- modelList$data
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 1)
  
  if(class(formula.) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(formula., "~~"))) else
    
    vars <- gsub(" ", "", unlist(strsplit(deparse(formula.), "~")))
  
  vars <- strsplit(vars, ":|\\*")
  
  if(!all(unlist(vars) %in% colnames(data)))
    
    stop("Variables not found in the model list. Ensure spelling is correct and remove all transformations!")
  
  residModList <- getResidModels(vars, modelList, data)
  
  if(all(class(residModList$ymod) == "numeric"))
    
    yresid <- data.frame(.id = names(residModList$ymod), yresid = residModList$ymod) else
      
      yresid <- data.frame(.id = rownames(GetData(residModList$ymod)), yresid = as.numeric(resid(residModList$ymod))) #resid.lme(ymod)
  
  if(all(class(residModList$xmod) == "numeric"))
    
    xresid <- data.frame(.id = names(residModList$xmod), xresid = residModList$xmod) else
      
      xresid <- data.frame(.id = rownames(GetData(residModList$xmod)), xresid = as.numeric(resid(residModList$xmod))) #resid.lme(xmod)
  
  rdata <- merge(yresid, xresid, by = ".id", all = TRUE)
  
  rdata <- rdata[order(as.numeric(as.character(rdata$.id))), -1]
  
  rownames(rdata) <- NULL
  
  return(rdata)
  
}

#' Calculate partial correlations from partial residuals
#' 
#' @keywords internal
#' 
partialCorr <- function(formula., modelList, data = NULL) {
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  if(is.null(data) & class(modelList) == "psem") data <- modelList$data
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 1)
  
  rdata <- partialResid(formula., modelList, data)
  
  rcor <- cor(rdata[, 1], rdata[, 2], use = "complete.obs")
  
  if(class(formula.) == "formula.cerror") vars <- gsub(" " , "", unlist(strsplit(formula., "~~"))) else
    
    vars <- gsub(" ", "", unlist(strsplit(deparse(formula.), "~")))
  
  vars <- strsplit(vars, ":|\\*")
  
  flag <- unlist(vars) %in% unlist(sapply(listFormula(modelList), function(x) all.vars.merMod(x)[1]))
  
  if(all(flag == FALSE)) {
    
    ctest <- cor.test(rdata[, 1], rdata[, 2])
    
    t. <- ctest$statistic
    
    N <- ctest$parameter
    
    P <- ctest$p.value
    
  } else {
    
    N <- nrow(rdata)
    
    residModList <- getResidModels(vars, modelList, data)
    
    k <- sum(sapply(residModList, function(x)
      
      if(all(class(x) == "numeric")) 0 else
        
        length(all.vars.merMod(formula(x))) - 2
      
    ) )
    
    k <- k[!duplicated(k)]
    
    k <- k[!k %in% vars]
    
    k <- length(k)
    
    N <- N - k - 2
    
    t. <- rcor * sqrt(N/(1 - rcor^2))
    
    P <- 1 - pt(t., N)
    
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
#' 
#' @keywords internal
#' 
resid.lme <- function(model) {
  
  if(any(class(model) %in% c("lme", "glmmPQL"))) {
    
    Q <- length(summary(model)$modelStruct$reStruct)
    
    r <- resid(model, level = 0:Q)
    
    r <- r[, 1]
    
  } else r <- resid(model)
  
  return(r)
  
}

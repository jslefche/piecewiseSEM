get.sem.fit = function(modelList, data, corr.errors = NULL, add.vars = NULL, 
                       grouping.vars = NULL, top.level.vars = NULL, adjust.p = FALSE, 
                       basis.set = NULL, pvalues.df = NULL, disp.conditional = FALSE,
                       .progressBar = TRUE) {

  if(!all(sapply(modelList, function(i) 
    all(class(i) %in% c("lm", "glm", "negbin", "lme", "lmerMod", "merModLmerTest", "glmerMod", "glmmPQL")) ) ) )
    stop("Model classes in model list are not supported")
  
  if(is.null(data)) stop("Must supply dataset to function")
  
  if(!all(unlist(lapply(modelList, nobs)))) 
    warning("All models do not have the same number of observations")
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, corr.errors, add.vars)
    
    basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 
    
  }
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")

  if(is.null(pvalues.df)) {
    
    if(.progressBar == T & length(basis.set) > 1) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
    
    pvalues.df = pvalues.df = get.missing.paths(modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, adjust.p, basis.set, disp.conditional, .progressBar) }
  
  fisher.c = get.fisher.c(modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, adjust.p, basis.set, pvalues.df, disp.conditional, .progressBar)
    
  AIC.c = get.aic(modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, adjust.p, basis.set, pvalues.df, disp.conditional, .progressBar)
  
  l = list(pvalues.df, fisher.c, AIC.c)
  
  names(l) = c("missing.paths", "Fisher.C", "AIC")

  return(l)
  
} 
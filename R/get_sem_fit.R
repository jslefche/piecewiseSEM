get.sem.fit = function(modelList, add.vars = NULL, adjust.p = FALSE, corr.errors = NULL, .progressBar = FALSE) {

  if(!all(sapply(modelList, function(i) 
    all(class(i) %in% c("lm", "glm", "negbin", "lme", "lmerMod", "merModLmerTest", "glmerMod", "glmmPQL")) ) ) )
    stop("Model classes in model list are not supported")
  
  if(!all(unlist(lapply(modelList, nobs)))) 
    warning("All models do not have the same number of observations")
  
  basis.set = get.basis.set(modelList, add.vars, corr.errors)
  
  basis.set = filter.exogenous(modelList, basis.set, add.vars)
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")

  pvalues.df = get.missing.paths(modelList, adjust.p, .progressBar, basis.set, add.vars, corr.errors)
  
  fisher.c = get.fisher.c(modelList, pvalues.df, adjust.p, .progressBar, basis.set, add.vars, corr.errors)
    
  AIC.c = get.aic(modelList, pvalues.df, add.vars, corr.errors, adjust.p, .progressBar, basis.set)
  
  l = list(pvalues.df, fisher.c, AIC.c)
  
  names(l) = c("missing.paths", "Fisher.C", "AIC")

  return(l)
  
} 
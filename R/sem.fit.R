sem.fit = function(
  
  modelList, data, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, top.level.vars = NULL, 
  adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL, .progressBar = TRUE
  
  ) {

  if(!all(sapply(modelList, function(i) 
    
    all(class(i) %in% c("lm", "glm", "negbin", "gls", "pgls", "lme", "lmerMod", "merModLmerTest", "glmerMod", "glmmPQL")) 
    
    ) ) ) stop("Model classes in model list are not supported")
  
  if(is.null(data)) stop("Must supply dataset")
  
  if(!all(unlist(lapply(modelList, nobs)))) warning("All models do not have the same number of observations")
  
  # Get basis set
  if(is.null(basis.set)) basis.set = sem.basis.set(modelList, corr.errors, add.vars)
  
  # Filter exogenous variables
  basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars)
  
  # Conduct d-sep tests
  if(is.null(pvalues.df)) pvalues.df = sem.missing.paths(
    
    modelList, data, conditional = FALSE, corr.errors, add.vars, grouping.vars,
    top.level.vars, adjust.p, basis.set, model.control, .progressBar
    
  )
  
  # Derive Fisher's C statistic and compare to Chi-squared distribution
  fisher.c = sem.fisher.c(
    
    modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, 
    adjust.p, basis.set, pvalues.df, model.control, .progressBar
    
  )
  
  # Use Fisher's C to derive AIC values
  AIC.c = sem.aic(
    
    modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, 
    adjust.p, basis.set, pvalues.df, model.control, .progressBar
    
  )
  
  # Round values in output table
  pvalues.df[, c(2:3, 5:6)] = apply(pvalues.df[, c(2:3, 5:6)], 2, function(x) round(x, 4) )
  
  # Return d-sep tests, Fisher's C, and AIC values in a list
  l = list(pvalues.df, fisher.c, AIC.c)
  
  names(l) = c("missing.paths", "Fisher.C", "AIC")

  return(l)
  
} 
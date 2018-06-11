sem.fit = function(
  
  modelList, data, conditional = FALSE, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, 
  grouping.fun = mean, adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL, 
  .progressBar = TRUE
  
  ) {
  
  if(is.null(data)) stop("Must supply dataset")
  
  n.obs = sapply(modelList, function(x) {
    
    if(class(x) == "rq") 
      
      length(na.omit(residuals(x))) else
        
        nobs(x)
    
  } )
  
  if(!all(n.obs)) warning("All models do not have the same number of observations")
  
  # Get basis set
  if(is.null(basis.set)) basis.set = suppressMessages(suppressWarnings(
    
    sem.basis.set(modelList, corr.errors, add.vars)
    
  ) ) 
  
  # Conduct d-sep tests
  if(is.null(pvalues.df)) pvalues.df = sem.missing.paths(
    
    modelList, data, conditional, corr.errors, add.vars, grouping.vars,
    grouping.fun, adjust.p, basis.set, model.control, .progressBar
  
  )
  
  # Derive Fisher's C statistic and compare to Chi-squared distribution
  fisher.c = sem.fisher.c(
    
    modelList, data, corr.errors, add.vars, grouping.vars, 
    grouping.fun, adjust.p, basis.set, pvalues.df, model.control, .progressBar
    
  )
  
  # Use Fisher's C to derive AIC values
  AIC.c = sem.aic(
    
    modelList, data, corr.errors, add.vars, grouping.vars,
    grouping.fun, adjust.p, basis.set, pvalues.df, model.control, .progressBar
    
  )
  
  # Round values in output table
  pvalues.df[, c(2:3, 5:6)] = apply(pvalues.df[, c(2:3, 5:6)], 2, function(x) round(x, 4) )
  
  # Return d-sep tests, Fisher's C, and AIC values in a list
  l = list(pvalues.df, fisher.c, AIC.c)
  
  names(l) = c("missing.paths", "Fisher.C", "AIC")

  return(l)
  
} 
sem.aic = function(
  
  modelList, data, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, top.level.vars = NULL, 
  adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL, .progressBar = TRUE
  
  ) {
  
  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))
  
  basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 
  
  if(is.null(pvalues.df)) pvalues.df = sem.missing.paths(
    
    modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars,
    adjust.p, basis.set, model.control, .progressBar
    
  )
  
  fisher.c = sem.fisher.c(
    
    modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, 
    adjust.p, basis.set, pvalues.df, model.control, .progressBar
    
  )
    
  # Calculate likelihood degrees of freedom  
  K = do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))
  
  # Calculate AIC
  AIC = unname(fisher.c[1] + 2 * K)
  
  # Calculate AICc
  AICc = unname(fisher.c[1] + 2 * K * (min(unlist(lapply(modelList, nobs)))/(min(unlist(lapply(modelList, nobs))) - K - 1)))
  
  # Return output in a data.frame
  data.frame(
    AIC = round(AIC, 3),
    AICc = round(AICc, 3), 
    K = round(K, 1), 
    n = round(min(unlist(lapply(modelList, nobs))), 1) )

}
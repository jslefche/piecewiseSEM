get.aic = function(modelList, corr.errors = NULL, add.vars = NULL, adjust.p = FALSE, 
                   basis.set = NULL, pvalues.df = NULL, .progressBar = TRUE) {
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, corr.errors, add.vars)
    
    basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 
    
  }
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")
  
  if(is.null(pvalues.df))
    pvalues.df = get.missing.paths(modelList, corr.errors, add.vars, adjust.p, basis.set, .progressBar)
  
  fisher.c = get.fisher.c(modelList, corr.errors, add.vars, adjust.p, basis.set, pvalues.df, .progressBar)
  
  K = do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))
  
  AIC = unname(fisher.c[1] + 2 * K)
  
  AICc = unname(fisher.c[1] + 2 * K * (max(unlist(lapply(modelList, nobs)))/(max(unlist(lapply(modelList, nobs))) - K - 1)))
  
  c(AIC = AIC, AICc = AICc, df = K)

}
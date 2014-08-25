get.aic = function(modelList, pvalues.df = NULL, add.vars = NULL, corr.errors = NULL, adjust.p = FALSE, 
                   .progressBar = FALSE, basis.set = NULL) {
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, add.vars, corr.errors)
    
    basis.set = filter.exogenous(modelList, basis.set, add.vars) 
    
  }
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")
  
  if(is.null(pvalues.df))
    pvalues.df = get.missing.paths(modelList, adjust.p, .progressBar, basis.set, add.vars, corr.errors)
  
  fisher.c = get.fisher.c(modelList, pvalues.df, adjust.p, .progressBar, basis.set, add.vars, corr.errors)
  
  K = do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))
  
  AIC = unname(fisher.c[1] + 2 * K)
  
  AICc = unname(fisher.c[1] + 2 * K * (max(unlist(lapply(modelList, nobs)))/(max(unlist(lapply(modelList, nobs))) - K - 1)))
  
  c(AIC = AIC, AICc = AICc, df = K)

}
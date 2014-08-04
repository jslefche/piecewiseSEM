get.aic = function(modelList, pvalues.df = NULL, adjust.p = FALSE, .progressBar = FALSE, basis.set = NULL, add.vars = NULL) {
  
  if(is.null(basis.set)) { 
    
    dag = dag.updated(modelList)
    
    basis.set = basiSet(dag)
    
    basis.set = lapply(basis.set, function(i) gsub(paste(LETTERS[1:10], collapse = ""), "\\:", i))
    
    basis.set = filter.exogenous(basis.set, add.vars) 
    
  }
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")
  
  if(is.null(pvalues.df)) {
    
    if(.progressBar == T & length(basis.set) > 1) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
    
    pvalues.df = get.missing.paths(modelList, adjust.p, .progressBar, basis.set, add.vars) }
  
  fisher.c = get.fisher.c(modelList, pvalues.df, adjust.p, .progressBar, basis.set, add.vars)
  
  K = do.call(sum, lapply(modelList, function(i) attr(logLik(i), "df")))
  
  AIC = fisher.c[1] + 2 * K
  
  AICc = fisher.c[1] + 2 * K * (mean(unlist(lapply(modelList, nobs)))/(mean(unlist(lapply(modelList, nobs))) - K - 1))
  
  c(AIC = AIC, AICc = AICc)

}
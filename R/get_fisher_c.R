get.fisher.c = function(modelList, pvalues.df = NULL, adjust.p = FALSE, .progressBar = FALSE, basis.set = NULL, add.vars = NULL) {
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, add.vars)
    
    basis.set = filter.exogenous(modelList, basis.set, add.vars) 
    
  }
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")

  if(is.null(pvalues.df)) {
    
    if(.progressBar == T & length(basis.set) > 1) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
    
    pvalues.df = get.missing.paths(modelList, adjust.p, .progressBar, basis.set, add.vars) }

  fisherC = -2 * sum(log(pvalues.df$p.value + 2e-16))
  
  p.value = if(is.null(pvalues.df)) NA else 1 - pchisq(fisherC, 2 * length(basis.set)) 
  
  c(Fisher.C = fisherC, P = p.value)
  
}
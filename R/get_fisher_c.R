get.fisher.c = function(modelList, data, corr.errors = NULL, add.vars = NULL, 
                        grouping.vars = NULL, top.level.vars = NULL, adjust.p = FALSE, 
                        basis.set = NULL, pvalues.df = NULL, disp.conditional = FALSE,
                        .progressBar = TRUE) {
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, corr.errors, add.vars)
    
    basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 
    
  }
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")

  if(is.null(pvalues.df)) {
    
    if(.progressBar == T & length(basis.set) > 1) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
    
    pvalues.df = pvalues.df = get.missing.paths(modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars,
                                                adjust.p, basis.set, disp.conditional, .progressBar) }

  fisherC = -2 * sum(log(pvalues.df$p.value + 2e-16))
  
  p.value = if(is.null(pvalues.df)) NA else 1 - pchisq(fisherC, 2 * length(basis.set)) 
  
  c(Fisher.C = fisherC, P = p.value)
  
}
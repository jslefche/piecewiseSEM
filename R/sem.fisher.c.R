sem.fisher.c = function(
  
  modelList, data, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, top.level.vars = NULL, 
  adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL, .progressBar = TRUE
  
  ) {
  
  if(is.null(basis.set)) basis.set = suppressWarnings(sem.basis.set(modelList, corr.errors, add.vars))
    
  basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 

  if(is.null(pvalues.df)) pvalues.df = suppressWarnings(sem.missing.paths(
    
    modelList, data, conditional = FALSE, corr.errors, add.vars, grouping.vars, 
    top.level.vars, adjust.p, basis.set, model.control, .progressBar
    
    ) )
  
  # Convert any p-values to a very small number as log(0) == -Inf
  if(length(basis.set) > 0 & any(pvalues.df$p.value == 0)) pvalues.df[pvalues.df$p.value == 0, "p.value"] = 2e-16
  
  # Calculate Fisher's C statistic
  fisher.C = if(length(basis.set) > 0) -2 * sum(log(pvalues.df$p.value)) else 0
  # Calculate associated p-value from Chi-squared distribution
  p.value = 1 - pchisq(fisher.C, 2 * length(basis.set)) 
  
  # Return output in a data.frame
  data.frame(fisher.c = round(fisher.C, 2), k = round(2 * length(basis.set), 1), p.value = round(p.value, 3))
  
}
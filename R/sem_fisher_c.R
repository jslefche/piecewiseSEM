sem.fisher.c = function(
  
  modelList, data, corr.errors = NULL, add.vars = NULL, grouping.vars = NULL, top.level.vars = NULL, 
  filter.exog = TRUE, adjust.p = FALSE, basis.set = NULL, pvalues.df = NULL, model.control = NULL, 
  .progressBar = TRUE
  
  ) {
  
  if(is.null(basis.set)) basis.set = sem.basis.set(modelList, corr.errors, add.vars)
    
  if(filter.exog == T) basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars) 

  if(is.null(pvalues.df)) pvalues.df = sem.missing.paths(
    
    modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, filter.exog,
    adjust.p, basis.set, model.control, .progressBar
    
    )
  
  # Convert any p-values to a very small number as log(0) == -Inf
  if(any(pvalues.df$p.value == 0)) pvalues.df[pvalues.df$p.value == 0, "p.value"] = 2e-16
  
  # Calculate Fisher's C statistic
  fisherC = if(length(pvalues.df) > 0) -2 * sum(log(pvalues.df$p.value)) else 0
  # Calculate associated p-value from Chi-squared distribution
  p.value = 1 - pchisq(fisherC, 2 * length(basis.set)) 
  
  # Return output in a data.frame
  data.frame(Fisher.C = round(fisherC, 2), k = round(length(basis.set), 1), P = round(p.value, 3))
  
}
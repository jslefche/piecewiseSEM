get.ppa.fit = function(modelList, data, corr.errors = NULL, add.vars = NULL, 
                       grouping.vars = NULL, top.level.vars = NULL, adjust.p = FALSE, 
                       basis.set = NULL, pvalues.df = NULL, disp.conditional = FALSE,
                       model.control = NULL, sig = 3, .progressBar = TRUE) {

  if(!all(sapply(modelList, function(i) 
    all(class(i) %in% c("lm", "glm", "negbin", "lme", "lmerMod", "merModLmerTest", "glmerMod", "glmmPQL","pgls")) ) ) )
    stop("Model classes in model list are not supported")
  
  if(is.null(data)) stop("Must supply dataset to function")
  
  if(!all(unlist(lapply(modelList, nobs)))) 
    warning("All models do not have the same number of observations")
  
  if(is.null(basis.set)) { 
    
    basis.set = get.basis.set(modelList, corr.errors, add.vars)
    basis.set = basis.set
    
  }
  
  if(is.null(pvalues.df)) {
    
    if(.progressBar == T & length(basis.set) > 1) pb = txtProgressBar(min = 0, max = length(basis.set), style = 3) else pb = NULL
    
    pvalues.df = get.condind.ppa(modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, adjust.p, basis.set, disp.conditional, model.control, .progressBar) }
  
  fisher.c = get.fisher.c(modelList, data, sig, corr.errors, add.vars, grouping.vars, top.level.vars, adjust.p, basis.set, pvalues.df, disp.conditional, model.control, .progressBar)
    
  AIC.c = get.aic(modelList, data, sig, corr.errors, add.vars, grouping.vars, top.level.vars, adjust.p, basis.set, pvalues.df, disp.conditional, model.control, .progressBar)
  
  pvalues.df[,c("estimate","std.error","crit.value","p.value")] = round(pvalues.df[,c("estimate","std.error","crit.value","p.value")], sig)
  
  l = list(pvalues.df, fisher.c, AIC.c)
  
  names(l) = c("missing.paths", "Fisher.C", "AIC")

  return(l)
  
} 
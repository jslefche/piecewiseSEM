shipley.test = function(modelList, add.vars = NULL, adjust.p = FALSE, .progressBar = FALSE) {

  if(!all(sapply(modelList, function(i) 
    all(class(i) %in% c("lm", "glm", "negbin", "lme", "lmerMod", "merModLmerTest", "glmerMod", "glmmPQL")) ) ) )
    stop("Model classes in model list are not supported")
  
  if(!all(unlist(lapply(modelList, nobs)))) 
    warning("All models do not have the same number of observations")
  
  dag = dag.updated(modelList, add.vars)
  
  basis.set = basiSet(dag)
  
  basis.set = lapply(basis.set, function(i) gsub(paste(LETTERS[1:10], collapse = ""), "\\:", i))
    
  basis.set = filter.exogenous(basis.set, add.vars)
  
  if(length(basis.set) < 1) 
    warning("All endogenous variables are conditionally dependent: no test of d-sep necessary")

  pvalues.df = get.missing.paths(modelList, adjust.p, .progressBar, basis.set, add.vars)
  
  fisher.c = get.fisher.c(modelList, pvalues.df, adjust.p, .progressBar, basis.set, add.vars)
    
  AIC.c = get.aic(modelList, pvalues.df, adjust.p, .progressBar, basis.set, add.vars)
    
  l = list(pvalues.df, fisher.c, AIC.c)
  
  names(l) = c("missing.paths", "Fisher.C", "AIC")

  return(l)
  
} 
get.partial.resid = function(y, x, modelList) {
  
  y = gsub(" ", "", y)
  x = gsub(" ", "", x)
  
  y.model = modelList[[match(y, unlist(lapply(modelList, function(i) as.character(formula(i)[2]))))]]
  
  if(is.null(y.model)) stop("Check spelling of correlated variables")
  
  if(all(x %in% as.character(formula(y.model)[[3]]))) stop("Y is a direct function of X, no partial residuals obtainable")
  
  y.nox.formula = gsub(x, "", paste(format(formula(y.model)), collapse = ""))
  
  if(all(class(y.model) %in% c("lme", "glmmPQL"))) 
    y.nox.model = update(y.model, fixed = formula(y.nox.formula)) else
      y.nox.model = update(y.model, formula = formula(y.nox.formula))
                   
  x.noy.formula = gsub(y, x, y.nox.formula)
  
  if(all(class(y.model) %in% c("lme", "glmmPQL"))) 
    x.noy.model = update(y.model, fixed = formula(x.noy.formula)) else
      x.noy.model = update(y.model, formula = formula(x.noy.formula))
  
  resids.data = data.frame(resid(y.nox.model), resid(x.noy.model) )

  names(resids.data)=c(paste(y, "given.others", sep="."), paste(x, "given.others", sep="."))
  
  return(resids.data)
  
}
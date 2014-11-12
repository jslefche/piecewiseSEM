get.partial.resid = function(formula = y ~ x, modelList, data, model.control = NULL) {
  
  vars = unlist(strsplit(deparse(formula), "~"))
  y = gsub(" ", "", vars[1])
  x = gsub(" ", "", vars[2])
  
  y.model = modelList[[match(y, unlist(lapply(modelList, function(i) as.character(formula(i)[2]))))]]
  
  if(is.null(y.model)) stop("Check spelling of correlated variables - must match exactly response in model formula!")
  
  if(all(x %in% strsplit(deparse(formula(y.model)[[3]]), "\\+")[[1]])) stop("Y is a direct function of X, no partial residuals obtainable")
  
  if(is.null(model.control)) {
    if(class(y.model) %in% c("lme", "glmmPQL")) control = lmeControl() else 
      if(class(y.model) %in% c("lmerMod", "merModLmerTest")) control = lmerControl() else
        if(class(y.model) %in% c("glmerMod")) control = glmerControl()} else {
          if(!is.null(model.control) & class(y.model) %in% c("lme", "glmmPQL"))
            control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), is.null))]] else
              if(!is.null(model.control) & class(y.model) %in% c("lmerMod", "merModLmerTest"))
                control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), function(x) any(x %in% "lmerControl")))]] else
                  if(!is.null(model.control) & class(y.model) %in% c("glmerMod"))
                    control = model.control[[which(sapply(lapply(model.control, function(x) attr(x, "class")), function(x) any(x %in% "glmerControl")))]] }
  
  if(all(class(y.model) %in% c("lme", "glmmPQL"))) {
    y.nox.formula = gsub(x, "", paste(deparse(y.model$call), collapse = ""))
    fixed = gsub(".*fixed = (.*), data = .*$", "\\1", y.nox.formula)
    random = gsub(".*random = (.*))", "\\1", y.nox.formula)
    y.nox.model = update(y.model, fixed = as.formula(fixed), random = as.formula(random), control = control, data = data) } else {
      y.nox.formula = gsub(x, "", paste(formula(y.model), collapse = ""))
      y.nox.model = update(y.model, formula = formula(y.nox.formula), control = control, data = data) }
  
  y.replace = gsub("\\(", "\\\\(", y)
  y.replace = gsub("\\)", "\\\\)", y.replace)
  
  x.noy.formula = gsub(y.replace, x, y.nox.formula)
  
  if(all(class(y.model) %in% c("lme", "glmmPQL"))) {
    fixed = gsub(".*fixed = (.*), data = .*$", "\\1", x.noy.formula)
    x.noy.model = update(y.model, fixed = formula(fixed), random = as.formula(random), control = control, data = data) } else 
      x.noy.model = update(y.model, formula = formula(x.noy.formula), control = control, data = data)
  
  resids.data = data.frame(resid(y.nox.model), resid(x.noy.model) )
  
  names(resids.data)=c(paste(y, "given.others", sep="."), paste(x, "given.others", sep="."))
  
  return(resids.data)
  
}
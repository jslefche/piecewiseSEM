get.partial.resid = function(.formula = y ~ x, modelList, model.control = NULL) {
  
  vars = unlist(strsplit(deparse(.formula), "~"))
  y = gsub(" ", "", vars[1])
  x = gsub(" ", "", vars[2])
  
  y.model = modelList[[match(y, unlist(lapply(modelList, function(i) gsub(".*\\((.*?)\\)+.*", "\\1", formula(i)[2]))))]]
  
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
  
  
  ### NEED TO FIX RANDOM FORMULA FOR BOTH NLME AND LME4
  
  random = as.formula(gsub(x, "", gsub(".*random = (.*))", "\\1",  paste(deparse(y.model$call), collapse = ""))))
  
  if(all(class(y.model) %in% c("lme", "glmmPQL")))
    y.nox.model = update(y.model, fixed = drop.terms(y.model$terms, grep(x, attr(y.model$terms, "term.labels")), keep.response = TRUE), random = random, control = control) else
      y.nox.model = update(y.model, formula = drop.terms(y.model$terms, grep(x, attr(y.model$terms, "term.labels")), keep.response = TRUE), control = control)
  
  x.sub = attr(y.model$terms, "term.labels")[grep(x, attr(y.model$terms, "term.labels"))][1]
  
  if(all(class(y.nox.model) %in% c("lme", "glmmPQL"))) 
    x.noy.model = update(y.model, fixed = reformulate(deparse(formula(y.nox.model)[[3]]), response = x.sub), random = random, control = control) else 
      x.noy.model = update(y.model, formula = reformulate(deparse(formula(y.nox.model)[[3]]), response = x.sub), control = control)
  
  ### 
  
  resids.data = data.frame(resid(y.nox.model), resid(x.noy.model) )
  
  names(resids.data)=c(paste(y, "given.others", sep="."), paste(x, "given.others", sep="."))
  
  return(resids.data)
  
}
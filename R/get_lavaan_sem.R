get.lavaan.sem = function(modelList, add.vars = NULL) {

  sem.model = paste(unlist(lapply(modelList, function(model) paste(format(formula(model)), collapse = ""))), collapse = "\n")

  if(!is.null(add.vars)) sem.model = paste(sem.model, 
    paste(unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))), collapse="\n"), sep = "\n")
  
  model = modelList[[which.max(lapply(modelList, nobs))]]
  
  if(all(class(model) %in% c("lm", "glm", "negbin"))) sem.model.data = model$sem.model else 
    if(all(class(model) %in% c("lme", "glmmPQL"))) sem.model.data = model$data else
      if(all(class(model) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) sem.model.data = model@frame
  
  sem(sem.model, sem.model.data) 

}
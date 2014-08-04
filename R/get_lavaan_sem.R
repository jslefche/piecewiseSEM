get.lavaan.sem = function(modelList, add.vars = NULL) {

  model = paste(unlist(lapply(modelList, function(i) paste(format(formula(i)), collapse = ""))), collapse = "\n")

  if(!is.null(add.vars)) model = paste(model, 
    paste(unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))), collapse="\n"), sep = "\n")
  
  if(all(class(i) %in% c("lm", "glm", "negbin"))) model.data = i$model else 
    if(all(class(i) %in% c("lme", "glmmPQL"))) model.data = i$data else
      if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) model.data = i@frame
  
  sem(model, model.data) 

}
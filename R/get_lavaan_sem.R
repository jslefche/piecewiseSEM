get.lavaan.sem = function(modelList, add.vars = NULL) {

  model = do.call(paste,c(lapply(modelList, function(i) 
    if(class(i) == "lme") 
      paste(formula(i)[[2]], formula(i)[[1]], formula(i)[[3]]) else
        paste(formula(i)[[2]], "~", paste(names(fixef(i))[-1], collapse = "+")) ), 
    sep = "\n"))
  
  if(!is.null(add.vars)) model = paste(model, 
    paste(unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))), collapse="\n"), sep = "\n")
  
  sem(model,data=if(class(modelList[[1]]) == "lme") modelList[[1]]$data else modelList[[1]]@frame) }
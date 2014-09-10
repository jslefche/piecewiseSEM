get.lavaan.sem = function(modelList, data, corr.errors = NULL, add.vars = NULL) {

  sem.model = paste(c(unlist(lapply(modelList, function(i) {
    
    if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) 
      form = format(nobars(formula(i))) else
        form = format(formula(i))
    
    paste(form, collapse = "") } ) ), corr.errors), collapse = "\n")

  if(!is.null(add.vars)) 
    
    sem.model = paste(sem.model, paste(unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))), collapse="\n"), sep = "\n")
  
  sem(sem.model, data)

}
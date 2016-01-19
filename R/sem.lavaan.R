sem.lavaan = function(modelList, data, compute.int = TRUE, corr.errors = NULL, add.vars = NULL, ...) {

  if(is.null(data)) stop("Must supply dataset to function")
  
  # Get list of formula from model list
  formula.list = lapply(modelList, function(i) 
    
    if(all(class(i) %in% c("lm", "glm", "negbin", "lme", "glmmPQL", "gls", "pgls"))) formula(i) else 
      
      if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) nobars(formula(i))
    
  )
  
  # Hand-calculate interactions
  if(compute.int == TRUE & any(grepl("\\*|\\:", unlist(formula.list)))) {
    
    formula.list = lapply(formula.list, function(i) {
      
      # Get rhs
      rhs = attr(terms(i), "term.labels")
      
      # Get interactions
      ints = rhs[grepl("\\:|\\*", rhs)]
      
      # Break apart interaction
      ints.split = sapply(ints, strsplit, "\\:")
      
      # Cycle through and add interactions to dataset
      if(length(ints) > 0) {
        
        for(j in ints.split) {
          
          data[, paste0(j, collapse = "_")] <<- with(data, eval(parse(text = paste(j, collapse = " * "))))
          
        } 
        
      }
        
      # Cycle through and remove interactions in formula list
      formula(paste0(all.vars(i)[1], " ~ ", paste0(gsub("\\:", "_", rhs), collapse = " + ")))
      
    } )
    
  }
      
  # Convert model formula to lavaan syntax
  sem.model = paste(formula.list, collapse = "\n")
  
  if(!is.null(corr.errors)) sem.model = paste(sem.model, paste(corr.errors, collapse = "\n"), collapse = "\n")

  if(!is.null(add.vars))
    
    sem.model = 
    paste(sem.model, 
          paste0(sapply(add.vars, function(x) as.formula(paste(x, x, sep = "~"))), collapse = "\n"),
          collapse = "\n")

  # Run lavaan SEM
  model = sem(sem.model, data, ...)

  return(model)
  
}
get.scaled.data = function(modelList, data, standardize) {

  # Identify variables that are transformed in the model formula
  transform.vars = unlist(lapply(modelList, function(i) {
    
    # Break apart formula
    vars = rownames(attr(terms(i), "factors"))
    
    # Identify transformed vars
    vars[grepl("log\\(|log10\\(|sqrt\\(|I\\(", vars)]
  
  } ) )
  
  # For each variables in transform.vars, perform transformation and store as original variable
  for(i in transform.vars) {
    
    # Get column name to transform
    col.nm = gsub("(.*)\\+.*", "\\1", gsub(".*\\((.*)\\)", "\\1", i))
    
    # Get column number
    col.no = which(colnames(data) == gsub(" ", "", col.nm))
    
    # Get actual transformation
    trsf = gsub("(.*)\\(.*\\)", "\\1", i)
    
    # Perform transformation
    data[, col.no] = eval(parse(text = gsub(col.nm, paste0("data[, ", col.no, "]"), i)))
    
  }
      
  # Get variables to scale, ignoring variables that are modeled to non-normal distributions
  vars.to.scale = unlist(lapply(modelList, function(i) {
    
    err = try(family(i), TRUE)
    
    if(grepl("Error", err[1]) | grepl("gaussian", err[1]))
      
      all.vars(formula(i)) else {
        
        warning(
          paste("Reponse '", formula(i)[2], "' is not modeled to a gaussian distribution: keeping response on original scale")
        )
        
        NULL }
    
  } ) )
  
  # Remove variables that are factors
  vars.to.scale = vars.to.scale[!vars.to.scale %in% colnames(data)[sapply(data, is.factor)]]
  
  # Remove duplicated variables
  vars.to.scale = vars.to.scale[!duplicated(vars.to.scale)]
  
  # Scale those variables by mean and SD, or by range
  if(class(data) == "comparative.data")
    
    newdata = data$data else newdata = data
  
  newdata[, vars.to.scale] = apply(newdata[, vars.to.scale], 2, function(x) 
    
    if(standardize == "scale") scale(x) else
      
      if(standardize == "range") (x-min(x, na.rm = T)) / diff(range(x, na.rm = T))
    
  )
  
  if(class(data) == "comparative.data") {
    
    data$data = newdata
    
    newdata = data
    
  }
  
  return(newdata)
  
  }
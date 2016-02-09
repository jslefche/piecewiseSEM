get.scaled.data = function(modelList, data, standardize) {

  if(any(sapply(modelList, class) == "pgls") | class(data) == "comparative.data") {
    
    # Extract data.frame
    newdata = data$data
  
  } else {
    
    newdata = data
    
  }
  
  # Identify variables that are transformed in the model formula
  transform.vars = unlist(lapply(modelList, function(i) {
    
    # Break apart formula
    if(class(i) == "pgls") vars = i$varNames else
      
      vars = rownames(attr(terms(i), "factors"))
    
    # Identify transformed vars
    vars[grepl("log\\(|log10\\(|sqrt\\(|I\\(", vars)]
  
  } ) )
  
  # Remove duplicates
  transform.vars = transform.vars[!duplicated(transform.vars)]
  
  # For each variables in transform.vars, perform transformation and store as original variable
  for(i in transform.vars) {
    
    # Get column name to transform
    col.nm = gsub("(.*)\\+.*", "\\1", gsub(".*\\((.*)\\)", "\\1", i))
    
    # Get column number
    col.no = which(colnames(newdata) == gsub(" ", "", col.nm))
    
    # Get actual transformation
    trsf = gsub("(.*)\\(.*\\)", "\\1", i)
    
    # Perform transformation
    newdata[, col.no] = eval(parse(text = gsub(col.nm, paste0("newdata[, ", col.no, "]"), i)))
    
  }
      
  # Get variables to scale, ignoring variables that are modeled to non-normal distributions
  vars.to.scale = unlist(lapply(modelList, function(i) {
    
    err = try(family(i), TRUE)
    
    if(grepl("Error", err[1]) | grepl("gaussian", err[1]))
      
      all.vars(formula(i)) else {
        
        warning(
          paste("Reponse '", formula(i)[2], "' is not modeled to a gaussian distribution: keeping response on original scale")
        )
        
        all.vars(formula(i))[-1] }
    
  } ) )
  
  # Remove variables that are factors
  vars.to.scale = vars.to.scale[!vars.to.scale %in% colnames(data)[sapply(data, function(x) any(is.factor(x) | is.character(x)))] ]
  
  # Remove duplicateds variables
  vars.to.scale = vars.to.scale[!duplicated(vars.to.scale)]
  
  # Scale those variables by mean and SD, or by range
  newdata[, vars.to.scale] = apply(newdata[, vars.to.scale], 2, function(x) {
    
    if(standardize == "scale") scale(x) else
      
      if(standardize == "range") (x-min(x, na.rm = T)) / diff(range(x, na.rm = T)) else
        
        x
    
  } )
  
  if(class(data) == "comparative.data") {
    
    data$data = newdata
    
    newdata = data
    
  }
  
  return(newdata)
  
  }
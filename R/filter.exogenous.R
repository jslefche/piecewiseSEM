filter.exogenous = function(modelList, basis.set, corr.errors = NULL, add.vars = NULL) {
  
  # Covnert model list into list of vectors
  modelFormulaList = lapply(modelList, function(i) {
    
    i = formula(i)
    
    if(grepl("cbind\\(.*\\)", paste(i[2]))) 
      
      v = c(paste(formula(i)[2]), all.vars(i)[-(1:2)]) else
        
        v = all.vars(i)
    
    return(v)
    
  } )
  
  # Get vector of predictor variables
  vars.cleanup = function(x) { x = x[!duplicated(x)]; x = gsub(" ", "" , x); return(x) }
  
  pred.vars = sapply(modelFormulaList, function(x) x[-1])
  
  pred.vars = vars.cleanup(pred.vars)
  
  # Get vector of response variables
  response.vars = sapply(modelFormulaList, function(x) x[1])
  
  response.vars = vars.cleanup(response.vars)
  
  # Get vector of variables that appear only as predictors and never as responses
  filter.vars = pred.vars[!pred.vars %in% response.vars]
  
  # Remove filtered variables when they appear as responses in the basis set
  basis.set = basis.set[!sapply(basis.set, function(i) any(i[2] %in% filter.vars))]
  
  # Ensure that no entry in the basis set is just the reverse of an existing claim
  basis.set. = lapply(basis.set, function(i) i[order(i)])
  
  basis.set = basis.set[!duplicated(basis.set.)]
  
  # Ensure that no entry in the basis set already exists in the model list
  basis.set = lapply(basis.set, function(i) 
    
    if(any(sapply(modelFormulaList, function(j) any(i[1:2] %in% j[1]) & all(i[1:2] %in% j)))) NULL else i
    
    )
  
  # Replace interaction : with * when interaction is response in basis set & remove entries attempting to predict interaction
  basis.set = lapply(basis.set, function(i) 
    
    if(is.null(i)) i else { if(i[2] %in% response.vars) gsub(":", "\\*", i) else NULL }
    
    )
  
  # Remove NULLs from basis set
  basis.set = basis.set[!sapply(basis.set, is.null)]
  
  if(length(basis.set) < 1) stop("All endogenous variables are conditionally dependent.\nTest of directed separation not possible!", call. = FALSE)
  
  return(basis.set)
  
}
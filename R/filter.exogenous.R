filter.exogenous = function(modelList, basis.set = NULL, corr.errors = NULL, add.vars = NULL) {
  
  if(is.null(basis.set)) basis.set = get.basis.set(get.sort.dag(get.formula.list(modelList, add.vars)))
  
  # Convert basis.set into list of formulae
  formulaList = lapply(basis.set, function(i) formula(paste0(i[2], " ~ ", paste0(i[-2], collapse = " + "))))

  # Get vector of predictor variables
  pred.vars = c(add.vars, unlist(lapply(modelList, function(x) all.vars(formula(x))[3])))
  
  pred.vars = pred.vars[!duplicated(pred.vars)]
  
  pred.vars = gsub(" ", "" , pred.vars)

  # Get vector of response variables
  response.vars = unlist(lapply(modelList, function(x) paste(formula(x)[2])))
  
  response.vars = response.vars[!duplicated(response.vars)]
  
  response.vars = gsub(" ", "" , response.vars)
  
  # Get vector of variables that appear only as predictors and never as responses
  filter.vars = pred.vars[!pred.vars %in% response.vars]
  
  # Remove filtered variables when they appear as responses in the basis set
  basis.set = basis.set[!sapply(basis.set, function(i) any(i[2] %in% filter.vars))]
  
  # Ensure that no entry in the basis set already exists in the model list
  modelFormulaList = lapply(modelList, function(i) {
    
    i = formula(i)
    
    if(grepl(",", i[2])) 
      
      c(paste0("cbind(", all.vars(i)[1], ",", all.vars(i)[2], ")"), all.vars(i)[-(1:2)]) else
        
        all.vars(i)
    
  } )
    
  basis.set = lapply(basis.set, function(i) 
    
    if(any(sapply(modelFormulaList, function(j) any(i[1:2] %in% j[1]) & all(i[1:2] %in% j)))) NULL else i
    
    )
  
  # Ensure that no entry in the basis set is just the reverse of an existing claim
  basis.set. = lapply(basis.set, function(i) i[order(i)])
  
  basis.set = basis.set[!duplicated(basis.set.)]
  
  # Replace interaction : with * when interaction is response in basis set & remove entries attempting to predict interaction
  basis.set = lapply(basis.set, function(i) 
    
    if(is.null(i)) i else { if(i[2] %in% response.vars) gsub(":", "\\*", i) else NULL }
    
    )
  
  # Remove NULLs from basis set
  basis.set = basis.set[!sapply(basis.set, is.null)]
  
  if(length(basis.set) < 1) warning("All endogenous variables are conditionally dependent.\nTest of directed separation not possible!", call. = FALSE)
  
  return(basis.set)
  
}
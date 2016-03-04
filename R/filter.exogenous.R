filter.exogenous = function(modelList, basis.set = NULL, corr.errors = NULL, add.vars = NULL) {
  
  # If not basis set, generate using sem.basis.set()
  if(is.null(basis.set)) basis.set = sem.basis.set(modelList, corr.errors, add.vars)

  # Convert basis.set into list of formulae
  formulaList = lapply(basis.set, function(i) formula(paste0(i[2], " ~ ", paste0(i[-2], collapse = " + "))))
  
  # Get vector of predictor variables
  pred.vars = c(add.vars, unlist(lapply(modelList, function(x) all.vars(formula(x))[-1])))
  
  pred.vars = pred.vars[!duplicated(pred.vars)]

  # Get vector of response variables
  response.vars = unlist(lapply(modelList, function(x) all.vars(formula(x))[1]))
  
  response.vars = response.vars[!duplicated(response.vars)]
  
  # Get vector of variables that appear only as predictors and never as responses
  filter.vars = pred.vars[!pred.vars %in% response.vars]
  
  # Remove filtered variables when they appear as responses in the basis set
  basis.set = basis.set[!sapply(formulaList, function(i) any(all.vars(i)[1] %in% filter.vars))]
    
  if(length(basis.set) < 1) warning("All endogenous variables are conditionally dependent.\nTest of directed separation not possible!", call. = FALSE)
  
  return(basis.set)
  
}
get.dag = function(formulaList) {
  
  # Insert placeholder for interaction symbol
  formulaList = lapply(formulaList, function(i) formula(gsub("\\:", "_____", Reduce(paste, deparse(i)))))

  # Strip transformations
  formulaList.new = lapply(formulaList, function(i) {
    
    v = all.vars(i)
    
    if(length(v) > 1) 
      
      formula(paste0(v[1], " ~ ", paste0(v[-1], collapse = " + "))) else
        
        formula(paste(v, " ~ ", v))
    
  } )
  
  # Get list of variables and also strip transformations
  vars = unlist(lapply(formulaList, function(i) {
    
    drop.vars = all.vars(i)[grepl("offset", rownames(attr(terms(i), "factors")))]
    
    all.vars(i)[!all.vars(i) %in% drop.vars]
    
  } ) )
  
  vars = unname(vars[!duplicated(vars)])
  
  # Create adjacency matrix
  amat = do.call(cbind, lapply(vars, function(i) {
 
    # Indentify formula where variable is response
    form.no = which(sapply(formulaList.new, function(j) all.vars(formula(j))[1] == i))
    
    if(length(form.no) == 0) rep(0, length(vars)) else {
      
      # Isolate variable from formula list
      form = formulaList.new[[which(sapply(formulaList.new, function(j) all.vars(formula(j))[1] == i))]]
      
      vars %in% all.vars(form) + 0
      
    }
    
  } ) )
  
  # Ensure diagonal is zero
  diag(amat) = 0
  
  # Name rows and columns
  dimnames(amat) = list(vars, vars)
  
  # Determine if graph is acylic
  if(!isAcyclic(amat)) warning("The graph contains directed cycles!")
  
  return(amat)
  
}
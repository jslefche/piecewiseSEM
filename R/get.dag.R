get.dag = function(formulaList) {
  
  # Expand interactions to include interaction and main effects
  formulaList = lapply(formulaList, function(i) 
    
    if(grepl("\\*|\\:", paste(format(formula(i)), collapse = ""))) {
      
      lhs = paste(rownames(attr(terms(i), "factors"))[1])
      
      rhs = attr(terms(i), "term.labels")
      
      # Collapse into formula
      rhs = paste(rhs, collapse = " + ")
      
      # Insert placeholder for interaction symbol :
      rhs = gsub("\\:", "_____", rhs)
      
      formula(paste(lhs, " ~ ", rhs))
      
    }
    
    else i
    
  )
  
  # Get list of variables
  vars = unique(unlist(lapply(formulaList, function(i)
    
    unlist(dimnames(attr(terms(i), "factors")))
  
    ) ) ) 
  
  # Create adjacency matrix
  amat = do.call(cbind, lapply(vars, function(i) {
 
    # Isolate variable from formula list
    form = formulaList[sapply(formulaList, function(j) all.vars(j)[1] == i)]
    
    if(length(form) == 0) rep(0, length(vars)) else 
      
      vars %in% attr(terms(form[[1]]), "term.labels") + 0
    
  } ) )
  
  # Ensure diagonal is zero
  diag(amat) = 0
  
  # Name rows and columns
  dimnames(amat) = list(vars, vars)
  
  # Determine if graph is acylic
  if(!isAcyclic(amat)) warning("The graph contains directed cycles!")
  
  return(amat)
  
}
get.dag = function(formulaList) {
  
  # Expand interactions to include interaction and main effects
  formulaList = lapply(formulaList, function(i) 
    
    if(grepl("\\*|\\:", paste(format(formula(i)), collapse = ""))) {
      
      lhs = paste(rownames(attr(terms(i), "factors"))[1])
      
      rhs = attr(terms(i), "term.labels")
      
      # Insert placeholder for interaction symbol :
      rhs = gsub("\\:", "_____", rhs)
      
      # Sort interactions so always alphabetical
      for(j in which(grepl("_____", rhs))) {
        
        # Split interactions and sort alphabetically
        int = unlist(lapply(strsplit(rhs[j], "_____"), sort))
          
        # Recombine 
        int.rec = paste(int, collapse = "_____")
          
        # Re-insert into formula
        rhs[j] = int.rec
          
        }
 
      # Collapse into formula
      rhs = paste(rhs, collapse = " + ")
      
      # And return full formula
      formula(paste(lhs, " ~ ", rhs))
      
    }
    
    else i
    
  )

  # Strip transformations
  formulaList.new = lapply(formulaList, function(i) {
    
    v = all.vars(i)
    
    formula(paste0(v[1], " ~ ", paste0(v[-1], collapse = " + ")))
    
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
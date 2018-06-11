get.dag = function(formulaList) {
  
  # Insert placeholder for interaction symbol
  formulaList = lapply(formulaList, function(i) formula(gsub("\\:", "_____", Reduce(paste, deparse(i)))))

  # Strip transformations
  formulaList.new = lapply(formulaList, function(i) {
    
    v = all.vars(i)
    
    if(grepl("cbind\\(.*\\)", paste(formula(i)[2]))) {
      
      v = c(paste(v[1], v[2], sep = ","), v[-(1:2)])
      
    }
    
    return(v)
    
  } )
  
  # Get list of variables and also strip transformations
  vars = unlist(lapply(formulaList, function(i) {
    
    drop.vars = all.vars(i)[grepl("offset", rownames(attr(terms(i), "factors")))]
    
    vars. = all.vars(i)[!all.vars(i) %in% drop.vars]
    
    # Check to see whether response is vector
    if(grepl("cbind\\(.*\\)", paste(formula(i)[2]))) {
      
      vars. = c(paste(vars.[1], vars.[2], sep = ","), vars.[-(1:2)])
      
    }
    
    return(vars.)
    
  } ) )
  
  vars = unname(vars[!duplicated(vars)])
  
  # Create adjacency matrix
  amat = do.call(cbind, lapply(vars, function(i) {
 
    # Identify formula where variable is response
    form.no = which(sapply(formulaList.new, function(j) j[1] == i))
    
    if(length(form.no) == 0) rep(0, length(vars)) else {
      
      # Isolate variable from formula list
      form = formulaList.new[[which(sapply(formulaList.new, function(j)j[1] == i))]]
      
      vars %in% form + 0
      
    }
    
  } ) )
  
  # Ensure diagonal is zero
  diag(amat) = 0
  
  # Name rows and columns
  dimnames(amat) = list(vars, vars)
  
  # Determine if graph is acylic
  if(all(colSums(amat) > 0)) stop("Model is non-recursive (cyclical)! Remove directed cycles and re-run.", call. = FALSE)
  
  return(amat)
  
}
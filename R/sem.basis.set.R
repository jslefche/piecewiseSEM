sem.basis.set = function(modelList, corr.errors = NULL, add.vars = NULL) {

  # Get list of formula from model list
  formulaList = get.formula.list(modelList, add.vars)
  
  # Stop if any response in the model list is the same as any other response
  if(any(duplicated(sapply(formulaList, function(x) all.vars(x)[1])))) 
    
    stop("Duplicate responses detected in the model list.\n
         Collapse multiple single regressions into a single multiple regression so that each response appears only once!")
  
  # Get adjacency matrix and sort by parent to child nodes
  amat = get.sort.dag(formulaList)
  
  # If intercept only model, add response variable to adjacency matrix
  if(any(unlist(lapply(modelList, function(i) deparse(formula(i)[2]) %in% c("~1", "~ 1"))))) {
    
    # Isolate intercept only model(s)
    responses = sapply(modelList[which(sapply(modelList, function(i) grepl("~ 1|~1", deparse(formula(i)))))],
                       
                       function(j) strsplit(paste(formula(j)), "~")[[2]]
                       
    )
    
    amat = cbind(
      
      rbind(amat, matrix(rep(0, dim(amat)[1]), nrow = 1, dimnames = list(responses))),
      
      matrix(rep(0, dim(amat)[1] + 1), ncol = 1, dimnames = list(NULL, responses))
      
    )
    
  }
  
  # Generate basis set
  if(all(amat == 0) & all(dim(amat) == 1)) basis.set = NULL else basis.set = get.basis.set(amat)
  
  # Replace placeholder for interaction symbol with :
  basis.set = lapply(basis.set, function(i) gsub(paste("_____", collapse = ""), "\\:", i))
  
  # Filter exogenous predictors from the basis set
  basis.set = filter.exogenous(modelList, basis.set, corr.errors, add.vars)

  # Re-apply transformations
  basis.set = lapply(basis.set, function(i) {
    
    # Get list of transformed predictors
    t.pvars = lapply(formulaList, function(x) colnames(attr(terms(x), "factors")))
    
    # Get list of untransformed predictors
    pvars = lapply(formulaList, function(i) {
      
      v = all.vars(i)[-1]
      
      if(grepl("c\\(.*\\)", paste(formula(i)[2]))) v = all.vars(i)[-(1:2)]
      
      return(v)
      
    } )
    
    # Re-transform predictors
    for(j in (1:length(i))[-2]) {
      
      # Get variable index for lookup
      idx = which(sapply(pvars, function(k) any(k %in% i[j])))
      
      idx. = unlist(sapply(pvars, function(l) which(l == i[j])))
      
      # Extract from transformed predictors
      i[j] = t.pvars[[idx]][idx.]
      
    }
    
    # Repeat for responses
    t.rvars = lapply(formulaList, function(x) rownames(attr(terms(x), "factors"))[1])
    
    # Get list of untransformed predictors
    rvars = lapply(formulaList, function(i) {
      
      v = all.vars(i)[1]
      
      if(grepl("c\\(.*\\)", paste(formula(i)[2]))) v = all.vars(i)[(1:2)]
      
      return(v)
      
    } )
    
    # Re-transform responses
    idx = which(sapply(rvars, function(k) any(k %in% i[2])))
    
    idx. = unlist(sapply(rvars, function(l) which(l == i[2])))
    
    # Extract from transformed responses
    i[2] = t.rvars[[idx]][idx.]
    
    return(i)
    
  } )
  
  # If correlated errors are present, remove them from the basis set
  if(!is.null(corr.errors)) {
    
    basis.set =  lapply(1:length(basis.set), function(i) {
      
      inset = unlist(lapply(corr.errors, function(j) {
        
        corr.vars = gsub(" ", "", unlist(strsplit(j,"~~")))
        
        all(
          
          unlist(
            
            lapply(1:2, function(k)
              
              grepl(paste(corr.vars, collapse = "|"), basis.set[[i]][k]) 
              
            ) 
          ) 
        ) 
        
      } ) )
      
      if(any(inset == TRUE)) NULL else basis.set[[i]]  
      
    } )
  }
  
  # Replace any d-sep where interactions are regressed against the main effect with NULL
  basis.set = lapply(basis.set, function(i) {
    
    if(is.null(i)) NULL else {
      
      if(grepl("\\*", i[1])) {
        
        int = strsplit(i[1], "\\*")[[1]]
        
        if(any(int %in% i[2])) NULL else i 
        
      } 
      
      else i
      
    }
    
  } )
  
  # Identify responses for which offset is present
  rpl = do.call(rbind, lapply(formulaList, function(i) {
    
    lhs = paste(rownames(attr(terms(i), "factors"))[1])
    
    rhs = rownames(attr(terms(i), "factors"))[-1]
    
    if(any(grepl("offset", rhs)))
      
      data.frame(response = lhs, offset = rhs[grepl("offset", rhs)]) else
        
        NULL
    
  } ) )
  
  # Add offset to basis set
  if(!is.null(rpl)) 
    
    basis.set = lapply(basis.set, function(i) {
      
      if(any(i[2] == rpl$response)) {
        
        c(i, as.character(rpl[rpl$response == i[2], "offset"]))
        
      } else i 
      
    } )
  
  # Remove NULLs from basis set
  basis.set = basis.set[!sapply(basis.set, is.null)]

  return(basis.set)
  
}
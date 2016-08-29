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

  
  ### TEMPORARY FIX ###
  
  # Reverse intermediate endogenous variables fitted to non-normal distributions
  basis.set = endogenous.reverse(basis.set, modelList)
  
  ### TEMPORARY FIX ###
  
    
  # Re-apply transformations
  basis.set = lapply(basis.set, function(i) {
    
    # Get list of transformed predictors
    t.pvars = lapply(formulaList, function(x) colnames(attr(terms(x), "factors")))
    
    # Get list of untransformed predictors
    pvars = lapply(formulaList, function(i) {
      
      if(grepl("cbind\\(.*\\)", paste(formula(i)[2]))) 
        
        v = all.vars(i)[-(1:2)] else
          
          v = all.vars(i)[-1]
        
        return(v)
        
    } )
    
    # Get list of transformed responses
    t.rvars = lapply(formulaList, function(x) gsub(" " , "", rownames(attr(terms(x), "factors"))[1]))
    
    # Get list of untransformed responses
    rvars = lapply(formulaList, function(i) {
      
      if(grepl("cbind\\(.*\\)", paste(formula(i)[2]))) 
        
        v = paste0("cbind(", paste(all.vars(i)[1:2], collapse = ","), ")") else
          
          v = all.vars(i)[1]
        
        return(gsub(" " , "", v))
        
    } )
    
    # Re-transform predictors
    for(j in (1:length(i))[-2]) {
      
      # Get variable index for lookup
      idx = cbind(
        
        which(sapply(pvars, function(k) any(k %in% i[j]))),
        
        unlist(sapply(pvars, function(l) which(l == i[j])))
        
      )
      
      if(sum(idx) > 0) {
        
        # Conduct lookup
        t.pvar = sapply(1:nrow(idx), function(m) t.pvars[[idx[m, 1]]][idx[m, 2]] )
        
        t.pvar = t.pvar[!duplicated(t.pvar)]
        
        # Replace predictors
        if(length(t.pvar) > 0) i[j] = t.pvar[which.max(sapply(t.pvar, function(p) nchar(p)))]
        
      }
      
    }
    
    # Get variable index for lookup
    idx = cbind(
      
      which(sapply(rvars, function(k) any(k %in% i[2]))),
      
      unlist(sapply(rvars, function(l) which(l == i[2])))
      
    )
    
    if(sum(idx) != 0) {
      
      # Conduct lookup
      t.rvar = sapply(1:nrow(idx), function(m) t.rvars[[idx[m, 1]]][idx[m, 2]] )
      
      t.rvar = t.rvar[!duplicated(t.rvar)]
      
      # Replace predictors
      if(length(t.rvar) > 0) i[2] = t.rvar[which.max(sapply(t.rvar, function(p) nchar(p)))]
      
    }
    
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
  
  # Re-assign names from dropped entries
  names(basis.set) = as.numeric(as.factor(names(basis.set)))
  
  return(basis.set)
  
}
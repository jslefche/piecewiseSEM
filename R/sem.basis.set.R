sem.basis.set = function(modelList, corr.errors = NULL, add.vars = NULL) {
  
  # Get list of formula from model list
  formulaList = lapply(modelList, function(i) 
    
    if(all(class(i) %in% c("lm", "glm", "negbin", "lme", "glmmPQL", "gls", "pgls", "glmmadmb"))) formula(i) else 
      
      if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) nobars(formula(i))
    
  )
  
  if(any(unlist(lapply(formulaList, is.null)))) stop("At least one model class not yet supported")
  
  # Check to see if any variables in the formula list appear in add.vars
  if(any(unlist(lapply(formulaList, all.vars)) %in% add.vars)) stop("Variables in the model list appear in add.vars!")
  
  # If additional variables are present, add them to the basis set
  if(!is.null(add.vars)) {
    
    # If interactions are specified with an asterisk, replace with semicolon
    add.vars = sapply(add.vars, function(x) gsub(" \\* ", "\\:", x))
    
    formulaList = append(formulaList, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep = "~")))))
    
  }
  
  # Generate adjacency matrix
  amat = get.dag(formulaList)
  
  # Sort adjacency matrix by parent to child nodes
  amat = get.sort.dag(amat)
  
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
    
    i.new = i
    
    # Find response in original formula list
    form = formulaList[[which(sapply(formulaList, function(x) all.vars(x)[1] == i[2]))]]
    
    # Get transformed vars for that response
    transform.vars = rownames(attr(terms(form), "factors"))
    
    # Get stripped vars for that response
    stripped.vars = all.vars(form)
    
    # Index over untransformed vars and replace
    i.new[which(i %in% stripped.vars)] = transform.vars[which(stripped.vars %in% i)]
    
    # Extract vector of predictors from all models
    preds.list = lapply(formulaList, all.vars)
    
    # Find formulae that contain untransformed variables
    for(j in which(!i.new %in% transform.vars)) {
      
      # Get list of transformed variables
      j.new = lapply(formulaList, function(k) {
        
        pred.vars = all.vars(k)[-1]
        
        # Sort interactions so always alphabetical
        for(l in which(grepl("\\:", pred.vars))) {
          
          # Split interactions and sort alphabetically
          int = unlist(lapply(strsplit(pred.vars[j], "\\:"), sort))
          
          # Recombine 
          int.rec = paste(int, collapse = ":")
          
          # Re-insert into formula
          pred.vars[l] = int.rec
          
        }
        
        colnames(attr(terms(k), "factors"))[pred.vars == i.new[j]]
        
      } )
      
      # Name and order transformed variables by adjacency matrix
      names(j.new) = sapply(formulaList, function(x) all.vars(x)[1])
      
      j.new = j.new[colnames(amat)]
      
      # Choose first instance where variable is not null
      new.var = unlist(j.new[sapply(j.new, function(x) length(x) > 0)])[1]
      
      if(!is.null(new.var)) i.new[j] = new.var
      
    }
    
    return(i.new)
    
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
      
      if(grepl("\\:", i[1])) {
        
        int = strsplit(i[1], "\\:")[[1]]
        
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
  
  if(length(basis.set) < 1) warning("All endogenous variables are conditionally dependent.\nTest of directed separation not possible!", call. = FALSE)
  
  return(basis.set)
  
}
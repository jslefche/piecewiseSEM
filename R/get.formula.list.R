get.formula.list = function(modelList, add.vars = NULL) {
  
  # Get list of formula from model list
  formulaList = lapply(modelList, function(i) 
    
    if(all(class(i) %in% c("lm", "rq", "glm", "negbin", "lme", "glmmPQL", "gls", "pgls", "glmmadmb"))) formula(i) else 
      
      if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod", "glmmTMB"))) nobars(formula(i))
    
  )
  
  if(any(unlist(lapply(formulaList, is.null)))) stop("At least one model class not yet supported", .call = FALSE)
  
  # Check to see if any variables in the formula list appear in add.vars
  if(any(unlist(lapply(formulaList, all.vars)) %in% add.vars)) stop("Variables in the model list appear in add.vars!")
  
  # If additional variables are present, add them to the basis set
  if(!is.null(add.vars)) {
    
    formulaList = append(formulaList, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep = "~")))))
    
  }
  
  # Expand interactions in formula list
  formulaList = lapply(formulaList, function(i) 
    
    if(grepl("\\*|:", paste(format(formula(i)), collapse = ""))) {
      
      lhs = paste(rownames(attr(terms(i), "factors"))[1])
      
      rhs = attr(terms(i), "term.labels")
      
      # Sort interactions so always alphabetical
      for(j in which(grepl(":", rhs))) {
        
        # Split interactions and sort alphabetically
        int = unlist(lapply(strsplit(rhs[j], ":"), sort))
        
        # Recombine 
        int.rec = paste(int, collapse = ":")
        
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
  
  if(any(sapply(formulaList, function(x) grepl("poly", x)))) 

      stop("Polynomials computed within the regression are not yet allowed.\nCompute externally and supply each component as a separate variable!")
  
  return(formulaList)
  
}
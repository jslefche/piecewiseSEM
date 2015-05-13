filter.exogenous = function(modelList, basis.set = NULL, corr.errors = NULL, add.vars = NULL) {
  
  # If not basis set, generate using sem.basis.set()
  if(is.null(basis.set)) basis.set = sem.basis.set(modelList, corr.errors, add.vars)
   
  # Get vector of predictor variables
  pred.vars = unique(
    
    c(add.vars,
      
      unlist(
      
        lapply(modelList, function(i) 
          
          if(all(class(i) %in% c("gls", "pgls"))) attr(coef(i), "names",) else
            
            attr(terms(i), "term.labels")
        )
      )
    )
  )
  
  # Get vector of response variables
  response.vars = unlist(
    
    lapply(modelList, function(i)
    
        if(all(class(i) %in% c("gls", "pgls"))) i$namey else
      
            rownames(attr(terms(i), "factors"))[1]
    )
  )
  
  # Get vector of variables that appear only as predictors and never as responses
  filter.vars = pred.vars[!pred.vars %in% c(response.vars, add.vars)]
  
  # Remove filtered variables when they appear as responses in the basis set
  basis.set = lapply(basis.set, function(i) 
    
    if(i[2] %in% filter.vars | any(i[1] %in% gsub(".*\\((.*)\\).*", "\\1", i[2:length(i)]))) NULL else i
    
  )
  
  basis.set[!sapply(basis.set, is.null)] 

}
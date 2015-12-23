sem.basis.set = function(modelList, corr.errors = NULL, add.vars = NULL) {
  
  # Get list of formula from model list
  formula.list = lapply(modelList, function(i) 
    
    if(all(class(i) %in% c("lm", "glm", "negbin", "lme", "glmmPQL", "gls", "pgls"))) formula(i) else 
    
        if(all(class(i) %in% c("lmerMod", "merModLmerTest", "glmerMod"))) nobars(formula(i))
    
    )
      
  if(any(unlist(lapply(formula.list, is.null)))) stop("At least one model class not yet supported")
  
  # If additional variables are present, add them to the basis set
  if(!is.null(add.vars)) 
    
    formula.list = append(formula.list, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep = "~")))))
  
  # Generate adjacency matrix
  amat = get.dag(formula.list)

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
  
  # Add offsets back in for given response
  # Identify responses for which offset is present
  rpl = do.call(rbind, lapply(formula.list, function(i) {
    
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
  
  # Replace edit in DAG() function in the ggm package
  # body(DAG)[[2]] = substitute(f <- list(...))
  
  if(length(basis.set) < 1) warning("All endogenous variables are conditionally dependent.\nTest of directed separation not possible!", call. = FALSE)
  
  return(basis.set)
  
}
sem.basis.set = function(modelList, corr.errors = NULL, add.vars = NULL) {
  
  # Get DAG from model list
  dag = lapply(modelList, function(i) 
    
    if(all(class(i) %in% c("lm", "glm", "negbin", "lme", "glmmPQL", "gls", "pgls"))) formula(i) else 
    
        if(all(class(i) %in% c("lmer", "glmerMod", "lme4"))) nobars(formula(i))
    
    )
      
  if(any(unlist(lapply(dag, is.null)))) stop("At least one model class not yet supported")
  
  # If additional variables are preesnt, add them to the basis set
  if(!is.null(add.vars)) dag = 
      
      append(dag, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))))
  
  # Expand interactions to include interaction and main effects
  dag = lapply(dag, function(i) 
    
    if(grepl("\\*|\\:", paste(format(formula(i)), collapse = ""))) {
      
      lhs = paste(rownames(attr(terms(i), "factors"))[1])
      
      rhs = paste(attr(terms(i), "term.labels"), collapse = " + ")
      
      # Insert placeholder for interaction symbol :
      rhs = gsub("\\:", "%%", rhs)
      
      formula(paste(lhs, " ~ ", rhs))
      
    }
    
    else i
    
  )
    
  # Modify DAG() function from ggm package to accept a list
  body(DAG)[[2]] = substitute(f <- dag) 

  # Generate adjacency matrix
  mat = DAG(dag)
  
  # If intercept only model, add response variable to adjacency matrix
  if(any(unlist(lapply(modelList, function(i) grepl("~ 1|~1", deparse(formula(i))))))) {
    
    # Isolate intercept only model(s)
    responses = sapply(modelList[which(sapply(modelList, function(i) grepl("~ 1|~1", deparse(formula(i)))))],
           
           function(j) strsplit(paste(formula(j)), "~")[[2]]
           
    )
    
    mat = cbind(
      
      rbind(mat, matrix(rep(0, dim(mat)[1]), nrow = 1, dimnames = list(responses))),
      
      matrix(rep(0, dim(mat)[1] + 1), ncol = 1, dimnames = list(NULL, responses))
      
    )

  }
  
  # Generate basis set
  if(all(mat == 0) & all(dim(mat) ==1 )) basis.set = NULL else basis.set = basiSet(mat)
  
  # Replace placeholder for interaction symbol with :
  basis.set = lapply(basis.set, function(i) gsub(paste("%%", collapse = ""), "\\:", i))
  
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
  
  # Remove NULLs from basis set
  basis.set = basis.set[!sapply(basis.set, is.null)]
  
  # Replace edit in DAG() function in the ggm package
  body(DAG)[[2]] = substitute(f <- list(...))
  
  if(length(basis.set) < 1) warning("All endogenous variables are conditionally dependent.\nTest of directed separation not possible!", call. = FALSE)
  
  return(basis.set)
  
}
get.basis.set = function(modelList, corr.errors = NULL, add.vars = NULL) {
  
  dag = lapply(modelList, function(i) 
    if(all(class(i) %in% c("lm", "glm", "negbin", "lme", "glmmPQL"))) formula(i) else 
      nobars(formula(i)) )
  
  if(is.null(add.vars)) 
    dag = dag else 
      dag = append(dag, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))))
  
  dag = lapply(dag, function(i) 
    if(grepl("\\*|\\:", paste(format(formula(i)), collapse = ""))) {
      f = paste(rownames(attr(terms(i), "factors"))[1], "~",paste(colnames(attr(terms(i), "factors")), collapse = "+"))
      f = gsub("\\:", paste("%*%", collapse = ""), f)
      formula(f) 
    } else i )
  
  body(DAG)[[2]] = substitute(f <- dag) 
  
  basis.set = basiSet(DAG(dag))
  
  if(length(basis.set) < 1) stop("All endogenous variables are conditionally dependent: model is satured.\n  Test of directed separation not possible!")
  
  basis.set = lapply(basis.set, function(i) gsub(paste(".\\%\\*\\%.", collapse = ""), "\\:", i))
  
  if(!is.null(corr.errors)) {
    
    basis.set =  lapply(1:length(basis.set), function(i) {
      
      inset = unlist(lapply(corr.errors, function(j) {
        
        corr.vars = gsub(" ", "", unlist(strsplit(j,"~~")))
        
        all(unlist(lapply(1:2, function(k)
          grepl(paste(corr.vars, collapse = "|"), basis.set[[i]][k]) ) ) ) 
        
        } ) )
        
      if(any(inset == TRUE)) NULL else basis.set[[i]]  
      
    } )
  }
  
  basis.set =  lapply(1:length(basis.set), function(i) {
    
    if(is.null(basis.set[[i]])) NULL else {
      
      if(grepl("\\:",basis.set[[i]][1])) {
        
        int = strsplit(basis.set[[i]][1],"\\:")[[1]]
        
        if(any(int %in% basis.set[[i]][2])) NULL else basis.set[[i]] } else basis.set[[i]]
      
    }
    
  } )
  
  basis.set = basis.set[!sapply(basis.set, is.null)]
  
  body(DAG)[[2]] = substitute(f <- list(...))
  
  return(basis.set)
  
}
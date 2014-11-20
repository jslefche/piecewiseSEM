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
  
  basis.set = lapply(basis.set, function(i) gsub(paste("\\%\\*\\%", collapse = ""), "\\:", i))
  
  if(!is.null(corr.errors)) {
  
    basis.set =  lapply(1:length(basis.set), function(i) {
      
      inset = unlist(lapply(corr.errors, function(j) {
          
        corr.vars = gsub(" ", "", unlist(strsplit(j,"~~")))
        
        basis.set.sub = c() 
        for(k in 1:2) basis.set.sub[k] = gsub(".*\\((.*?)\\)+.*", "\\1", basis.set[[i]][k])
      
        all(basis.set.sub %in% corr.vars) } ))
      
      if(any(inset == TRUE)) NULL else basis.set[[i]]  
        
      } )
    }
  
  body(DAG)[[2]] = substitute(f <- list(...))
  
  basis.set = lapply(basis.set, function(i) gsub(" ", "", i))
  
  basis.set = basis.set[!sapply(basis.set, is.null)]
  
  return(basis.set)
  
}
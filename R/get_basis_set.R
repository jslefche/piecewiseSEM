get.basis.set = function(modelList, add.vars = NULL, corr.errors = NULL) {
  
  dag = lapply(modelList, function(i) 
    if(all(class(i) %in% c("lm", "glm", "negbin", "lme", "glmmPQL"))) formula(i) else 
      nobars(formula(i)) )
  
  if(is.null(add.vars)) 
    dag = dag else 
      dag = append(dag, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))))
  
  dag = lapply(dag, function(i) if(grepl("\\*|\\:", paste(format(formula(i)), collapse = ""))) {
    f = paste(formula(i)[[2]], "~", paste(colnames(attr(terms(i), "factors")), collapse = "+"))
    f = gsub("\\:", paste(LETTERS[1:10], collapse = ""), f)
    formula(f) }  else i )
   
  body(DAG)[[2]] = substitute(f <- dag) 
  
  basis.set = basiSet(DAG(dag))
  
  basis.set = lapply(basis.set, function(i) gsub(paste(LETTERS[1:10], collapse = ""), "\\:", i))
  
  if(!is.null(corr.errors)) {
  
    basis.set =  lapply(1:length(basis.set), function(i) {
      
      inset = unlist(lapply(corr.errors, function(j) {
          
      corr.vars = gsub(" ", "", unlist(strsplit(j,"~~")))
      
      all(basis.set[[i]][1:2] %in% corr.vars) } ))
      
      if(any(inset == TRUE)) NULL else basis.set[[i]]  
        
      } )
    }
  
  body(DAG)[[2]] = substitute(f <- list(...))
  
  return( basis.set[!sapply(basis.set, is.null)] )
  
}
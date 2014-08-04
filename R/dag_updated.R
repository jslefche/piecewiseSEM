dag.updated = function(modelList, add.vars = NULL) {
  
  dag = lapply(modelList, function(i) 
    if(all(class(i) %in% c("lm", "glm", "negbin", "lme", "glmmPQL"))) formula(i) else 
      nobars(formula(i)) )
  
  if(is.null(add.vars)) 
    dag = dag else 
      dag = append(dag, unname(sapply(add.vars, function(x) as.formula(paste(x, x, sep="~")))))
  
  dag = lapply(dag, function(i) if(grepl("\\*", format(formula(i)))) {
    f = paste(formula(i)[[2]], "~", paste(colnames(attr(terms(i), "factors")), collapse = "+"))
    f = gsub("\\:", paste(LETTERS[1:10], collapse = ""), f)
    formula(f) }  else i )
   
  body(DAG)[[2]] = substitute(f <- dag) 
  
  dag = DAG(dag)

  body(DAG)[[2]] = substitute(f <- list(...))
  
  return(dag)
  
}
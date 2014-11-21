filter.exogenous = function(modelList, basis.set = NULL, corr.errors = NULL, add.vars = NULL) {
  
  if(is.null(basis.set)) basis.set = get.basis.set(modelList, corr.errors, add.vars)
  
  pred.vars = unique(c(add.vars, unlist(lapply(modelList, function(i) attr(terms(i), "term.labels")))))
  
  response.vars = unlist(lapply(modelList, function(i) rownames(attr(terms(i), "factors"))[1]))

  filter.vars = pred.vars[!pred.vars %in% response.vars]
  
  basis.set = lapply(1:length(basis.set), function(i) 
    if(basis.set[[i]][2] %in% filter.vars | 
#          all(basis.set[[i]] %in% pred.vars) | 
         any(basis.set[[i]][1] %in% gsub(".*\\((.*)\\).*", "\\1", basis.set[[i]][2:length(basis.set[[i]])]))) NULL else 
      basis.set[[i]])
  
  basis.set[!sapply(basis.set, is.null)] 

}
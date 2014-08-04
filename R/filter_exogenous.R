filter.exogenous = function(basis.set, add.vars = NULL) {
  
  exogenous.vars = c(add.vars, unlist(lapply(modelList, function(i) colnames(attr(terms(i), "factors")))) )
  
  pred.vars = unlist(lapply(modelList, function(i) rownames(attr(terms(i), "factors"))[1]))
  
  filter.vars = exogenous.vars[!exogenous.vars %in% pred.vars]
  
  basis.set = lapply(1:length(basis.set), function(i) if(basis.set[[i]][2] %in% filter.vars) NULL else basis.set[[i]])
  
  basis.set[!sapply(basis.set, is.null)] 

}
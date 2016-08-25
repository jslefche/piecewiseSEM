endogenous.reverse = function(basis.set, modelList) {
  
  # Identify d-sep tests among endogenous variables and reverse if family != "gaussian"
  names(basis.set) = 1:length(basis.set)
  
  # Get sorted adjacency matrix
  amat = get.sort.dag(get.formula.list(modelList))
  
  # Identify intermediate endogenous variables
  idx = colnames(amat)[which(colSums(amat[which(colSums(amat) == 0), , drop = FALSE]) > 0)]
  
  idx = idx[!idx %in% names(which(colSums(amat[which(!colSums(amat) == 0), , drop = FALSE]) > 0))]
  
  # Identify variables in the basis set where intermediate endogenous variables are the response
  idx. = sapply(modelList, function(x) all.vars(formula(x))[1] %in% idx)
  
  if(any(idx.)) {
    
    # Get index of models for which responses are in index
    idx.. = idx[sapply(modelList[idx.], function(x) any(class(x) %in% c("glm", "negbin", "glmmPQL", "glmerMod")))]
    
    if(length(idx..) > 0) {
      
      basis.set = 
        
        append(basis.set, 
               
               lapply(
                 
                 which(sapply(basis.set, function(i) i[2] %in% idx..)), 
                 
                 function(i) 
                   
                   c(basis.set[[i]][2], basis.set[[i]][1], basis.set[[i]][-(1:2)])
                 
               )
               
        )
        
    }
    
  }
  
  return(basis.set)
  
}
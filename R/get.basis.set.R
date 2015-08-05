get.basis.set = function(amat) {
  
  ret = lapply(1:ncol(amat), function(i) {
    
    lapply(1:nrow(amat), function(j) {
      
      # If any relationship is directional or the variable predicts itself, return NULL
      if(amat[j, i] != 0 | amat[i, j] != 0 | i == j) NULL else {
        
        # Get directed relationship
        dsep = unlist(dimnames(amat[j, i, drop = FALSE]))
        
        
        
        
        # Determine if direct relationship is cylic
        if(all(rowSums(amat[, dsep, drop = FALSE]) > 0)) NULL else {
        
          
          
          
          
          # Get vector of conditional variables
          cond.var = rownames(amat[, i, drop = FALSE])[which(amat[, i, drop= FALSE] == 1)]
          
          c(
            unlist(dimnames(amat[j, i, drop = FALSE])),
            cond.var[!cond.var %in% rownames(amat[j, , drop = FALSE])]
          )
          
        }
        
      }
      
    } ) 
    
  } )
  
  ret = unlist(ret, recursive = FALSE)
  
  ret = ret[!sapply(ret, is.null)]
  
  return(ret)
  
}
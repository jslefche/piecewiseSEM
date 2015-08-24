get.basis.set = function(amat) {
  
  # Sort adjacency matrix by parent to child nodes
  amat = topSort(amat)
  
  ret = lapply(1:ncol(amat), function(j) {
    
    lapply(j:nrow(amat), function(i) {
      
      if(amat[j, i] != 0 | j == i) NULL else {  
        
        # Get variables for independence test
        dsep = unlist(dimnames(amat[j, i, drop = FALSE]))
        
        # Get vector of conditional variables
        cond.var = c(
          rownames(amat)[which(amat[, dsep[1], drop = FALSE] == 1)],
          rownames(amat)[which(amat[, dsep[2], drop = FALSE] == 1)]
        )
        
        # Remove conditional variables already in the independence claim
        cond.var = cond.var[!cond.var %in% dsep]
        
        # Return full independence claim
        c(dsep, cond.var)
        
      }
      
    } ) 
    
  } )
  
  ret = unlist(ret, recursive = FALSE)
  
  ret = ret[!sapply(ret, is.null)]
  
  return(ret)
  
}
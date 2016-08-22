get.basis.set = function(amat) {
  
  ret = lapply(1:ncol(amat), function(i) {
    
    lapply(i:nrow(amat), function(j) {
      
      if(amat[i, j] != 0 | i == j) NULL else {  
        
        # Get variables for independence test
        dsep = unlist(dimnames(amat[i, j, drop = FALSE]))
        
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

  ret = lapply(ret, function(j) j[!duplicated(j)])
    
  ret = ret[!sapply(ret, is.null)]
  
  # Add binding for binomial variables
  for(i in 1:length(ret)) {
    
    if(any(grepl(",", ret[[i]]))) {
      
      idx = which(grepl(",", ret[[i]]))
      
      for(j in idx) {
        
        ret[[i]][j] = paste0("cbind(", ret[[i]][j], ")")
        
      }
      
    }
    
  }
  
  return(ret)
  
}
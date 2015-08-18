get.basis.set = function(amat) {
  
  ret = lapply(1:ncol(amat), function(i) {
    
    lapply(1:nrow(amat), function(j) {
      
      # If any relationship is directional or the variable predicts itself, return NULL
      if(amat[j, i] != 0 | amat[i, j] != 0 | i == j) NULL else {
        
        # Get directed relationship
        dsep = unlist(dimnames(amat[j, i, drop = FALSE]))
    
        # Determine if direct relationship is cylical
        resp = names(amat[, dsep[1]][amat[, dsep[1]] > 0])
        
        A = FALSE
        
        while(A == FALSE) {

          nm = names(amat[, resp[length(resp)]][amat[, resp[length(resp)]] > 0])
          
          if(length(nm) ==0) A = TRUE else resp = c(resp, nm)
          
        }
          
        # If relationship is cyclical then return NULL
        if(any(resp %in% dsep[-1])) NULL else {
          
          # Get vector of conditional variables
          cond.var = c(
            rownames(amat)[which(amat[, dsep[1], drop = FALSE] == 1)],
            rownames(amat)[which(amat[, dsep[2], drop = FALSE] == 1)]
          )
       
          # Return full independence claim
          c(
            dsep,
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
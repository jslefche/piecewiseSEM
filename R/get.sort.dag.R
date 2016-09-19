get.sort.dag = function(formulaList) {

  # Get adjaceny matrix
  amat = get.dag(formulaList)
  
  # Get predictors where colSums == 0
  idx = unname(which(colSums(amat) == 0))

  # Of remaining variables, look for only those with links to the predictors in col.zero
  idx. = 
    which(!colnames(amat) %in% colnames(amat)[idx])[
      colSums(amat[!colnames(amat) %in% colnames(amat)[idx], !colnames(amat) %in% colnames(amat)[idx], drop = FALSE]) == 0 & 
        colSums(amat[colnames(amat) %in% colnames(amat)[idx], !colnames(amat) %in% colnames(amat)[idx], drop = FALSE]) != 0]
  
  # Sort by increasing position in amat
  pos = apply(amat[, idx., drop = FALSE], 2, function(x) max(which(x == 1)))
    
  idx. = c(idx, which(colnames(amat) %in% names(pos[order(pos)])))
 
  # Sort remaining variables by increasing position in amat
  pos. = apply(amat[, !colnames(amat) %in% colnames(amat)[idx.], drop = FALSE], 2, function(x) max(which(x == 1)))
  
  idx.. = c(idx., which(colnames(amat) %in% names(pos.[order(pos.)])))

  # Return adjacency matrix, sorted
  amat[idx.., idx..]
  
}
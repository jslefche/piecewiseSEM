sort.dag = function(amat) {

  # Get predictors where colSums == 0
  col.zero = colnames(amat)[colSums(amat) == 0]

  # Put those predictors at the beginning of the index
  idx = col.zero
  
  # Of remaining variables, look for only those with links to the predictors in col.zero
  idx = c(idx, 
          colnames(amat[, !colnames(amat) %in% idx, drop = FALSE])[
            colSums(amat[!rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) == 0 & 
              colSums(amat[rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 ]
  )
  
  # Of remaining variables, look for those with links in both col.zero and !col.zero
  idx = c(idx, 
          colnames(amat[, !colnames(amat) %in% idx, drop = FALSE])[
            colSums(amat[!rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 & 
              colSums(amat[rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 ]
  )

  # Finally, only variables with links in !col.zero
  idx = c(idx, 
          colnames(amat[, !colnames(amat) %in% idx, drop = FALSE])[
            colSums(amat[!rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 & 
              colSums(amat[rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) == 0 ]
  )
             
  # Return adjacency matrix, sorted
  amat[idx, idx]
  
}
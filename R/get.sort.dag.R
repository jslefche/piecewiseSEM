get.sort.dag = function(formulaList) {

  # Get adjaceny matrix
  amat = get.dag(formulaList)
  
  # Get predictors where colSums == 0
  col.zero = colnames(amat)[colSums(amat) == 0]

  # Put those predictors at the beginning of the index
  idx = col.zero
  
  # Of remaining variables, look for only those with links to the predictors in col.zero
  first.pred = colnames(amat[, !colnames(amat) %in% idx, drop = FALSE])[
    colSums(amat[!rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) == 0 & 
      colSums(amat[rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 ]
  
  first.pred = names(sort(colSums(amat[, first.pred, drop = FALSE]), decreasing = TRUE))
  
  idx = c(idx, first.pred)
  
  # Of remaining variables, look for those with links in both col.zero and !col.zero
  second.pred = colnames(amat[, !colnames(amat) %in% idx, drop = FALSE])[
    colSums(amat[!rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 & 
      colSums(amat[rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 ]

  second.pred = names(sort(colSums(amat[, second.pred, drop = FALSE]), decreasing = TRUE))
  
  idx = c(idx, second.pred)

  # Finally, only variables with links in !col.zero
  final.pred = colnames(amat[, !colnames(amat) %in% idx, drop = FALSE])[
    colSums(amat[!rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) != 0 & 
      colSums(amat[rownames(amat) %in% col.zero, !colnames(amat) %in% idx, drop = FALSE]) == 0 ]

  final.pred = names(sort(colSums(amat[, final.pred, drop = FALSE]), decreasing = TRUE))
  
  idx = c(idx, final.pred)
             
  # Return adjacency matrix, sorted
  amat[idx, idx]
  
}
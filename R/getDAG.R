#' Generate adjacency matrix from list of structural equations
#'
#' @param modelList A list of structural equations
#' 
#' @export
#' 

getDAG <- function(modelList) {
  
  el <- getEdgelist(modelList)
  
  graph <- graph_from_edgelist(el)
  
  amat <- as.matrix(as_adjacency_matrix(graph))
  
  order <- names(topo_sort(graph))
  
  amat <- amat[order, order]
  
  if(
    sum(colSums(amat) > 0) < 1 &
    is_dag(graph)
  ) 
    
    stop("Model is non-recursive. Remove feedback loops!", call. = FALSE)
  
  
  return(amat)
  
}


#' Get edge list from list of equations
#' 
#' @keywords internal
#' 
getEdgelist <- function(modelList) {
  
  fList <- listFormula(modelList, formulas = 0)
  
  do.call(rbind, lapply(fList, function(i) {
    
    f <- all.vars_notrans(i)
    
    matrix(c(rep(f[1], length(f) - 1), f[-1]), ncol = 2)
    
  } ) )
  
}

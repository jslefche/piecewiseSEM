#' Generate adjacency matrix from list of structural equations
#'
#' @param modelList A list of structural equations
#'
#' @export
#'
getDAG <- function(modelList) {

  formulaList <- listFormula(modelList)

  fList <- lapply(formulaList, function(i) if(any(class(i) %in% "formula.cerror")) NULL else i)

  fList <- fList[!sapply(fList, is.null)]

  fList <- lapply(fList, all.vars_trans, smoothed = TRUE)

  vars <- unlist(fList)

  vars <- unname(vars[!duplicated(vars)])

  amat <- do.call(cbind, lapply(vars, function(i) {

    f <- which(sapply(fList, function(j) j[1] == i))

    if(length(f) == 0) rep(0, length(vars)) else {

      form <- fList[[which(sapply(fList, function(k) k[1] == i))]]

      vars %in% form + 0

    }

  } ) )

  dimnames(amat) <- list(vars, vars)

  diag(amat) <- 0

  amat <- sortDag(amat, fList)

  if(!igraph::is_dag(igraph::graph_from_adjacency_matrix(amat)))
    
    stop("Model is non-recursive. Remove feedback loops!", call. = FALSE)

  return(amat)

}

#' Sort DAG based on ancestry
#'
#' @keywords internal
#'
sortDag <- function(amat, formulaList) { 
  
  counter <- sapply(rownames(amat), function(i) {
    
    indicated <- i
    
    flag <- TRUE
    
    counter <- 0
    
    while(flag) {
      
      pos <- sapply(formulaList, function(k) k[1] %in% indicated)
      
      if(all(!pos) || sum(amat[, i]) == 0) flag <- FALSE 
      
      else { 
        
        counter <- counter + 1
        
        indicated <- formulaList[pos][[1]][-1]
        
        if(counter == 1) flag <- TRUE
        
        else{
          
          if(any(indicated %in% memory_indicated)) { flag <- FALSE ; counter <- counter - 1}
          
          else flag <- TRUE
          
        }
        
        if(counter == 1) memory_indicated <- indicated
        
        else memory_indicated <- c(memory_indicated, indicated)
        
      }
      
    }
    
    return(counter)
    
  } )
  
  amat <- amat[names(sort(counter)), names(sort(counter))]
  
  return(amat)
  
}

# #' Generate adjacency matrix from list of structural equations
# #'
# #' @param modelList A list of structural equations
# #' 
# #' @import igraph
# #' 
# #' @export
# #' 
# getDAG <- function(modelList) {
# 
#   el <- getEdgelist(modelList)
# 
#   graph <- graph_from_edgelist(el)
# 
#   amat <- as.matrix(as_adjacency_matrix(graph))
# 
#   order <- names(topo_sort(graph))
# 
#   amat <- amat[order, order]
# 
#   if(
#     sum(colSums(amat) > 0) < 1 &
#     is_dag(graph)
#   )
# 
#     stop("Model is non-recursive. Remove feedback loops!", call. = FALSE)
# 
# 
#   return(amat)
# 
# }
# 
# #' Get edge list from list of equations
# #'
# #' @keywords internal
# #'
# getEdgelist <- function(modelList) {
# 
#   fList <- listFormula(modelList, formulas = 0)
# 
#   do.call(rbind, lapply(fList, function(i) {
# 
#     f <- all.vars_trans(i)
# 
#     matrix(c(f[-1], rep(f[1], length(f) - 1)), ncol = 2)
# 
#   } ) )
# 
# }
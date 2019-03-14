#' Generate adjacency matrix from list of structural equations
#'
#' @param modelList A list of structural equations
#' 
#' @export
#' 
Dag <- function(modelList) {
  
  formulaList <- listFormula(modelList)

  fList <- lapply(formulaList, function(i) if(any(class(i) %in% "formula.cerror")) NULL else i)

  fList <- fList[!sapply(fList, is.null)]

  fList <- lapply(fList, all.vars_trans)

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

  if(
    sum(colSums(amat) > 0) < 1 &
    igraph::is_dag(igraph::graph_from_adjacency_matrix(amat))
    ) 
    
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

    while(flag == TRUE) {

      pos <- sapply(formulaList, function(k) k[1] %in% indicated)

      if(all(pos == FALSE) | sum(amat[, i]) == 0) flag <- FALSE else {

        counter <- counter + 1

        indicated <- formulaList[pos][[1]][-1]

        flag <- TRUE

      }

    }

    return(counter)

  } )

  amat = amat[names(sort(counter)), names(sort(counter))]

  return(amat)

}

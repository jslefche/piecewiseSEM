#' Generate adjacency matrix from list of structural equations
#'
#' @param formulaList a list of formulae corresponding to structural equations

Dag <- function(formulaList) {

  fList <- lapply(formulaList, function(i) if(any(class(i) %in% "formula.cerror")) NULL else i)

  fList <- fList[!sapply(fList, is.null)]

  fList <- lapply(fList, all.vars.merMod)

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

  if(cyclic(amat, fList)) stop("Model is recursive. Remove feedback loops and re-run model")

  return(amat)

}

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

cyclic <- function(amat, formulaList) {

  vars <- colnames(amat[, colSums(amat) > 0, drop = FALSE])

  cyc = sapply(vars, function(i) {

    indicated <- i

    indicated_sum <- c()

    flag <- TRUE

    counter <- FALSE

    while(flag == TRUE) {

      pos <- sapply(formulaList, function(k) k[1] %in% indicated)

      if(all(pos == FALSE)) flag <- FALSE else {

        indicated <- formulaList[pos][[1]][-1]

        indicated_sum <- c(indicated_sum, indicated)

        flag <- TRUE

      }

    }

    i %in% indicated_sum

  } )

  any(cyc == TRUE)

}

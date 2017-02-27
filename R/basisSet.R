#' Get basis set for test of directed separation

#' @param modelList a list of structural equations

basisSet <- function(modelList) {

  formulaList <- listFormula(modelList)

  amat <- Dag(formulaList)

  b <- lapply(1:nrow(amat), function(i) {

    lapply(i:ncol(amat), function(j) {

      if(amat[i, j] != 0 | i == j) NULL else {

        cvars <- c(
          rownames(amat)[amat[, rownames(amat)[i], drop = FALSE] == 1],
          rownames(amat)[amat[, rownames(amat)[j], drop = FALSE] == 1]
        )

        cvars <- cvars[!duplicated(cvars)]

        c(rownames(amat)[i], rownames(amat)[j], cvars)

      }

    } )

  } )

  b <- unlist(b, recursive = FALSE)

  b <- filterExogenous(b, amat)

  b <- b[!sapply(b, is.null)]

  b <- removeCerror(b, formulaList)

  # Reverse intermediate endogenous variables

  b <- replaceTrans(modelList, b, amat)

  return(b)

}

#' Filter relationships among exogenous variables from the basis set

filterExogenous <- function(b, amat) {

  exo <- colnames(amat[, colSums(amat) == 0, drop = FALSE])

  b <- lapply(b, function(i)

    if(all(i[1:2] %in% exo)) NULL else i

  )

  return(b)

}

#' Remove correlated errors from the basis set

removeCerror <- function(b, formulaList) {

  ceList <- lapply(formulaList, function(i) if(any(class(i) == "formula.cerror")) {

    strsplit(i, " ~~ ")[[1]]

    } else NULL)

  ceList <- ceList[!sapply(ceList, is.null)]

  if(length(ceList) > 0) {

    b <- lapply(b, function(i) {

      lapply(ceList, function(j) {

        if(all(i[1:2] %in% j)) NULL else i

      } )

    } )

    b <- unlist(b, recursive = FALSE)

    b <- b[!sapply(b, is.null)]

  }

  return(b)

}

#' Reverse non-linear intermediate endogenous variables

reverseEndo <- function(b, modelList) {




}

#' Replace transformations in the basis set by cycling through neighbors and applying
#' transformations in order of how variables are treated in the child nearest to current node

replaceTrans <- function(modelList, b, amat) {

  formulaList <- listFormula(modelList)

  trans <- lapply(formulaList, all.vars.trans)

  notrans <- lapply(formulaList, all.vars.notrans)

  b <- lapply(b, function(i) {

    f <- which(sapply(notrans, function(j) j[1] == i[2]))

    if(length(i) > 2)

      i[3:length(i)] <- trans[[f]][which(notrans[[f]] %in% i[3:length(i)])]

    flag <- TRUE

    res <- i[2]

    while(flag == TRUE) {

      f <- which(sapply(notrans, function(j) res == j[1]))

      if(sum(f) == 0) flag <- FALSE else {

        if(i[1] %in% notrans[[f]][-1]) {

          i[1] <- trans[[f]][which(notrans[[f]] == i[1])]

          flag <- FALSE

          } else {

            res <- colnames(amat)[(which(res == colnames(amat)) - 1)]

          }

      }

    }

    return(i)

  } )

  return(b)

}

all.vars.trans <- function(.formula) {

  if(class(.formula) == "formula")

    rownames(attr(terms(.formula), "factors")) else

      strsplit(.formula, " ~~ ")[[1]]

}

all.vars.notrans <- function(.formula) {

  if(class(.formula) == "formula")

    all.vars(.formula) else

      strsplit(.formula, " ~~ ")[[1]]

}

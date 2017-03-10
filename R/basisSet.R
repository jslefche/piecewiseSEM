#' Get basis set for test of directed separation

#' @param modelList a list of structural equations

basisSet <- function(modelList, direction = NULL) {

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

  b <- reverseNonLin(modelList, b, amat, formulaList)

  b <- replaceTrans(modelList, b, amat)

  if(!is.null(direction))

    b <- specifyDir(b, direction)

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

    unlist(strsplit(i, " ~~ "))

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

#' Get vector of transformed variables
all.vars.trans <- function(.formula) {

  if(class(.formula) == "formula") {

    n <- rownames(attr(terms(.formula), "factors"))

    if(any(grepl("\\|", n))) {

      idn <- which(grepl("\\|", n))

      f <- rownames(attr(terms(.formula), "factors"))[-idn]

      return(f)

    } else return(n)

  } else unlist(strsplit(.formula, " ~~ "))

}

#' Get vector of untransformed variables
all.vars.notrans <- function(.formula) {

  if(class(.formula) == "formula")

    all.vars.merMod(.formula) else

      unlist(strsplit(.formula, " ~~ "))

}

#' If intermediate endogenous variables are nonlinear, return both directions
reverseNonLin <- function(modelList, b, amat, formulaList) {

  modelList <- modelList[sapply(formulaList, function(x) !class(x) %in% c("formula.cerror"))]

  formulaList <- formulaList[sapply(formulaList, function(x) !class(x) %in% c("formula.cerror"))]

  names(b) <- 1:length(b)

  idx <- which(colSums(amat[colSums(amat) == 0, , drop = FALSE]) > 0)

  idx <- idx[!idx %in% which(colSums(amat[!colSums(amat) == 0, , drop = FALSE]) > 0)]

  idx <- names(idx)

  idm <- sapply(formulaList, function(i) all.vars.merMod(i)[1] %in% idx)

  idm <- idm[

    sapply(modelList[idm], function(x) {

    .family <- try(family(x), silent = TRUE)

    if(class(.family) == "try-error") FALSE else TRUE

  } ) ]

  if(length(idm) > 0) {

    if(any(sapply(modelList[idm], function(x) family(x)$family != "gaussian"))) {

      idf <- idx[sapply(modelList[idm], function(x) family(x)$family != "gaussian")]

      if(length(idf) > 0) {

        b <- append(b, lapply(b[sapply(b, function(x) x[2] %in% idf)], function(i)

          c(i[2], i[1], i[-(1:2)]) )

          )

      }

    }

  }

  r <- sapply(formulaList, function(x) all.vars.merMod(x)[1])

  b <- b[sapply(b, function(x) any(x[2] %in% r))]

  return(b)

}

#' Replace transformations in the basis set by cycling through neighbors and applying
#' transformations in order of how variables are treated in the child nearest to current node

replaceTrans <- function(modelList, b, amat) {

  formulaList <- listFormula(modelList)

  trans <- lapply(formulaList, all.vars.trans)

  notrans <- lapply(formulaList, all.vars.notrans)

  b <- lapply(b, function(i) {

    f <- which(sapply(notrans, function(j) j[1] == i[2]))

    flag <- TRUE

    res <- i[2]

    while(flag == TRUE) {

      f <- which(sapply(notrans, function(j) res == j[1]))

      if(sum(f) == 0) flag <- FALSE else {

        for(j in f) {

          if(i[1] %in% notrans[j][-1]) {

            i[1] <- trans[j][which(notrans[j] == i[1])]

            flag <- FALSE

          } else {

            res <- colnames(amat)[(which(res == colnames(amat)) - 1)]

          }

        }

      }

    }

    return(i)

  } )

  return(b)

}

#' Remove items from the basis set whose direction is a priori specified
specifyDir <- function(b, direction) {

  vars <- gsub(" ", "", unlist(strsplit(direction, "\\->|<\\-")))

  b[which(sapply(b, function(i) i[1] == vars[2] & i[2] == vars[1]))] <- NULL

  return(b)

}

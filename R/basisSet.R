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

  b <- b[!sapply(b, is.null)]

  if(length(b) > 0) {

    b <- filterExisting(b, formulaList)

    b <- filterExogenous(modelList, b, amat)

    b <- filterInteractions(b)

    b <- removeCerror(b, formulaList)

    b <- reverseAddVars(modelList, b, amat)

    b <- reverseNonLin(modelList, b, amat)

    if(!is.null(direction)) b <- specifyDir(b, direction)

    # b <- replaceTrans(modelList, b, amat)

  }

  return(b)

}

#' Remove existing paths from the basis set
filterExisting <- function(b, formulaList) {

  b <- lapply(b, function(i) {

    f <- formulaList[sapply(formulaList, function(x) all.vars.trans(x)[1] == i[1])]

    if(any(sapply(f, function(x) any(all.vars.trans(x)[-1] %in% i[2])))) NULL else i

  } )

  b <- b[!sapply(b, is.null)]

  return(b)

}

#' Filter relationships among exogenous variables from the basis set (ignoring add.vars)
filterExogenous <- function(modelList, b, amat) {

  formulaList <- listFormula(modelList, formulas = 3)

  exo <- colnames(amat[, colSums(amat) == 0, drop = FALSE])

  exo <- exo[!exo %in% sapply(formulaList, function(x) ifelse(x[[3]] == 1, deparse(x[[2]]), ""))]

  b <- lapply(b, function(i)

    if(all(i[1:2] %in% exo) | i[2] %in% exo) NULL else i

  )

  b <- b[!sapply(b, is.null)]

  return(b)

}

#' Filter interactions from the d-sep tests
filterInteractions <- function(b) {

  b <- lapply(b, function(i) if(any(grepl("\\:", i[1:2]))) NULL else i )

  b <- b[!sapply(b, is.null)]

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

      flag = any(sapply(ceList, function(j) all(i[1:2] %in% j)))

      if(flag == TRUE) NULL else i

    } )

    b <- b[!sapply(b, is.null)]

  }

  return(b)

}

#' Replace transformations in the basis set by cycling through neighbors and applying
#' transformations in order of how variables are treated in the child nearest to current node
# replaceTrans <- function(modelList, b, amat) {
#
#   if(length(b) > 0) {
#
#     modelList <- removeData(modelList, formulas = 1)
#
#     formulaList <- listFormula(modelList)
#
#     amat <- amat[, colSums(amat) > 0]
#
#     formulaList <- formulaList[sapply(formulaList, function(x) which(all.vars.notrans(x)[1] == colnames(amat)))]
#
#     trans <- lapply(formulaList, all.vars.trans)
#
#     notrans <- lapply(formulaList, all.vars.notrans)
#
#     b <- lapply(rev(b), function(i) {
#
#       i[2] <- trans[[which(sapply(notrans, function(x) x[1] == i[2]))]][1]
#
#       flag <- TRUE
#
#       while(flag == TRUE) {
#
#         for(j in (1:length(i))[-2]) {
#
#           for(k in length(formulaList):1) {
#
#             if(i[j] %in% notrans[[k]][-1]) {
#
#               i[j] <- trans[[k]][which(i[j] == notrans[[k]])]
#
#               flag <- FALSE
#
#             }
#
#           }
#
#         }
#
#       }
#
#       return(i)
#
#     } )
#
#     b <- rev(b)
#
#   }
#
#   return(b)
#
# }

#' Reverse added variables (e.g., y ~ 1)
reverseAddVars <- function(modelList, b, amat) {

  formulaList <- listFormula(modelList, formulas = 3)

  exo <- colnames(amat[, colSums(amat) == 0, drop = FALSE])

  exo <- exo[exo %in% sapply(formulaList, function(x) ifelse(x[[3]] == 1, deparse(x[[2]]), ""))]

  b <- lapply(b, function(i) if(i[2] %in% exo) c(i[2], i[1], i[-(1:2)]) else i)

  return(b)

}

#' If intermediate endogenous variables are nonlinear, return both directions
reverseNonLin <- function(modelList, b, amat) {

  if(length(b) > 0) {

    modelList <- removeData(modelList, formulas = 1)

    formulaList <- listFormula(modelList, formulas = 1)

    names(b) <- 1:length(b)

    idx <- which(colSums(amat[colSums(amat) == 0, , drop = FALSE]) > 0)

    idx <- idx[!idx %in% which(colSums(amat[!colSums(amat) == 0, , drop = FALSE]) > 0)]

    idx <- names(idx)

    idm <- sapply(formulaList, function(i) all.vars.trans(i)[1] %in% idx)

    idm <- idm[

      sapply(modelList[idm], function(x) {

        .family <- try(family(x), silent = TRUE)

        if(class(.family) == "try-error") FALSE else

          if(.family$family == "gaussian") FALSE else

            TRUE

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

    r <- sapply(formulaList, function(x) all.vars.trans(x)[1])

    b <- b[sapply(b, function(x) any(x[2] %in% r))]

  }

  return(b)

}

#' Remove items from the basis set whose direction is a priori specified
specifyDir <- function(b, direction) {

  vars <- gsub(" ", "", unlist(strsplit(direction, "\\->|<\\-")))

  b[which(sapply(b, function(i) i[1] == vars[2] & i[2] == vars[1]))] <- NULL

  return(b)

}

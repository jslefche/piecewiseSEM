#' Derivation of the basis set
#' 
#' Acquires the set of independence claims--or the 'basis set'--for use in
#' evaluating the goodness-of-fit for piecewise structural equation models.
#' 
#' This function returns a list of independence claims. Each
#' claim is a vector of the predictor of interest, followed by the response,
#' and, if present, any conditioning variables.
#' 
#' Relationships among exogenous variables are omitted from the basis set
#' because the directionality is unclear--e.g., does temperature 
#' cause latitude or does latitude cause temperature?--and the assumptions 
#' of the variables are not specified in the list of structural equations, 
#' so evaluating the relationship becomes challenging without further 
#' input from the user. This creates a circular scenario whereby the 
#' user specifies relationships among exogenous variables, raising the 
#' issue of whether they should be included as directed paths if they can 
#' be assigned directional relationships.
#' 
#' Paths can be omitted from the basis set by specifying them as correlated
#' errors using \code{\%~~\%} or by assigning a directionality using 
#' the argument \code{direction}, e.g. \code{direction = c("X <- Y")}. 
#' This can be done if post hoc examination of the d-sep tests reveals
#' nonsensical independence claims (e.g., arthropod abundance predicting
#'  photosynthetically-active radiation) that the user may wish to 
#'  exclude from evaluation.
#' 
#' @param modelList A list of structural equations
#' @param direction a vector of claims defining the specific directionality of any independence 
#' claim(s) 
#' @param interactions whether interactions should be included in independence claims. 
#' Default is FALSE
#' 
#' @return A \code{list} of independence claims.
#' 
#' @author Jon Lefcheck <LefcheckJ@@si.edu>
#' @seealso \code{\link{dSep}}
#' 
#' @references Shipley, Bill. "A new inferential test for path models based on directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.
#' 
#' @export
#' 
basisSet <- function(modelList, direction = NULL, interactions = FALSE) {

  amat <- getDAG(modelList)

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
    
    b <- filterExogenous(b, modelList, amat)
    
    b <- filterSmoothed(b, modelList)
    
    b <- filterExisting(b, modelList)
    
    b <- filterInteractions(b, interactions)

    b <- removeCerror(b, modelList)

    b <- reverseAddVars(b, modelList, amat)

    b <- reverseNonLin(b, modelList, amat)

    b <- fixCatDir(b, modelList)
    
    if(!is.null(direction)) b <- specifyDir(b, direction)

    # b <- replaceTrans(modelList, b, amat)

  }
  
  class(b) <- "basisSet"

  return(b)

}

#' Remove existing paths from the basis set
#' 
#' @keywords internal
#' 
filterExisting <- function(b, modelList) {
  
  formulaList <- listFormula(modelList)
  
  b <- lapply(b, function(i) {

    f <- formulaList[sapply(formulaList, function(x) all.vars_trans(x)[1] == i[1])]

    if(any(sapply(f, function(x) any(all.vars_trans(x)[-1] %in% i[2])))) NULL else i

  } )

  b <- b[!sapply(b, is.null)]

  return(b)

}

#' Filter relationships among exogenous variables from the basis set (ignoring add.vars)
#' 
#' @keywords internal
#' 
filterExogenous <- function(b, modelList, amat) {

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
#' 
#' @keywords internal
#' 
filterInteractions <- function(b, interactions = FALSE) {

  if(interactions == FALSE) {
    
    b <- lapply(b, function(i) if(any(grepl("\\:", i[1:2]))) NULL else i )
    
  } else {
    
    b <- lapply(b, function(i) {
      
      vars <- unlist(sapply(i, strsplit, ":"))
      
      if(any(vars[1] %in% vars[-1]) | grepl(":", vars[2])) NULL else i
      
    } )
    
  }

  b <- b[!sapply(b, is.null)]

  return(b)

}

#' First, remove claims where linear and non-linear terms appear in the same claim
#' 
#' @keywords internal
#' 
filterSmoothed <- function(b, modelList) {
  
  b <- lapply(b, function(i) {
    
    vars <- gsub("(.*)\\,.*", "\\1", gsub("s\\((.*)\\).*", "\\1", i)) 
    
    if(vars[1] == vars[2] | any(vars[1] %in% vars[-(1:2)])) NULL else i
    
  } )
  
  b <- b[!sapply(b, is.null)]
  
  return(b)
  
}

#' Remove correlated errors from the basis set
#' 
#' @keywords internal
#' 
removeCerror <- function(b, modelList) {
  
  formulaList <- listFormula(modelList)

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
#' 
#' @keywords internal
#' 
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
#     formulaList <- formulaList[sapply(formulaList, function(x) which(all.vars_notrans(x)[1] == colnames(amat)))]
#
#     trans <- lapply(formulaList, all.vars_trans)
#
#     notrans <- lapply(formulaList, all.vars_notrans)
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

#' Reverse added variables
#' 
#' @keywords internal
#' 
reverseAddVars <- function(b, modelList, amat) {

  formulaList <- listFormula(modelList, formulas = 3)

  exo <- colnames(amat[, colSums(amat) == 0, drop = FALSE])

  exo <- exo[exo %in% sapply(formulaList, function(x) ifelse(x[[3]] == 1, deparse(x[[2]]), ""))]

  b <- lapply(b, function(i) if(i[2] %in% exo) c(i[2], i[1], i[-(1:2)]) else i)

  return(b)

}

#' If intermediate endogenous variables are nonlinear, return both directions
#' 
#' @keywords internal
#' 
reverseNonLin <- function(b, modelList, amat) {

  if(length(b) > 0) {

    modelList <- removeData(modelList, formulas = 1)

    formulaList <- listFormula(modelList, formulas = 1)

    names(b) <- 1:length(b)

    idx <- which(colSums(amat[colSums(amat) == 0, , drop = FALSE]) > 0)

    idx <- idx[!idx %in% which(colSums(amat[!colSums(amat) == 0, , drop = FALSE]) > 0)]

    idx <- names(idx)

    idm <- sapply(formulaList, function(i) all.vars_trans(i)[1] %in% idx)

    idm <- idm[

      sapply(modelList[idm], function(x) {

        family. <- try(family(x), silent = TRUE)

        if(inherits(family., "try-error") | all(is.na(family.))) FALSE else

          if(family.$family == "gaussian") FALSE else

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

    r <- sapply(formulaList, function(x) all.vars_trans(x)[1])

    b <- b[sapply(b, function(x) any(x[2] %in% r))]

  }

  return(b)

}

#' Remove duplicate items from the basis set whose direction is not a priori specified
#' 
#' @keywords internal
#' 
specifyDir <- function(b, direction) {
  
  vars <- gsub(" ", "", unlist(strsplit(direction, "\\->|<\\-")))
  
  b[which(sapply(b, function(i) i[1] == vars[2] & i[2] == vars[1]))] <- NULL
  
  return(b)
  
  # rels <- lapply(direction, function(d) { trimws(strsplit(d, "\\->|<\\-")[[1]]) })
  # 
  # dirs <- sapply(direction, function(d) { gsub(".*(\\->|<\\-).*", "\\1", d) })
  # 
  # for(i in 1:length(rels)) {
  #   
  #   fix <- flipOne(rels[[i]], dirs[i], b)
  #   
  #   b[[fix[[2]]]] <- fix[[1]] 
  #   
  #   }
  # 
  # return(b)

}


flipOne <- function(rel, arrow, b) {
  
  b_idx <- which(sapply(b, function(i) i[1] %in% rel & i[2] %in% rel))
  
  cond <- b[b_idx]

  if(arrow == "<\\-") { 
    
    r1 <- rel[1]
    
    rel[1] <- rel[2]
    
    rel[2] <- r1
    
    } 

  if(cond[[1]][1] != rel[1]) {
    
    cond[[1]][1] <- rel[1]
    
    cond[[1]][2] <- rel[2]
    
    }
  

  return(list(cond[[1]], b_idx))
  
}

#' Flip independence claims so categorical variables are not the response
#' 
#' @keywords internal
#' 
fixCatDir <- function(b, modelList) {
  
  b <- lapply(b, function(i) {
    
    var <- i[2]
    
    var <- gsub(".*\\((.*)\\).*", "\\1", var)
    
    data <- as.data.frame(modelList$data)
    
    if(class(data[, var]) %in% c("factor", "character")) {
      
      var1 <- i[1]
      
      i[1] <- var
      
      i[2] <- var1
      
    } 
    
    return(i)
    
  } )
  
  return(b)
  
}

#' Print basis set
#' 
#' @method print basisSet
#' 
#' @param x a basis set
#' @param ... further arguments passed to or from other methods
#' 
#' @export
#' 
print.basisSet <- function(x, ...) { 
  
  if(length(x) == 0) print("No independence claims in basis set.") else {
  
  ret <- lapply(x, function(oneLine) {
    
    st <- paste(oneLine[1], "|", oneLine[2], sep = " ")
    
    if(length(oneLine) > 2) st <- paste(st, "(", paste(oneLine[3:length(oneLine)], collapse = ", "), ")")
    
    st
    
  } )
  
  names(ret) <- paste("Claim", 1:length(ret))
  
  print(ret)
  
  }
  
}

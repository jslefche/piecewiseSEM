#' Tests of directed separation
#' 
#' Evaluation of conditional independence claims to be used in determining the
#' goodness-of-fit for piecewise structural equation models.
#' 
#' In cases involving non-normally distributed responses in the independence
#' claims that are modeled using generalized linear models, the significance of
#' the independence claim is not reversible (e.g., the P-value of Y ~ X is not
#' the same as X ~ Y). This is due to the transformation of the response via
#' the link function. In extreme cases, this can bias the goodness-of-fit
#' tests. \code{summary.psem} will issue a warning when this case is present
#' and provide guidance for solutions.
#' 
#' One solution is to specify the directionality of the relationship using the
#' \code{direction} argument, e.g. \code{direction = c("X <- Y")}. Another is
#' to run both tests (Y ~ X, X ~ Y) and return the most conservative (i.e.,
#' lowest) P-value, which can be toggled using the \code{conserve = TRUE}
#' argument.
#' 
#' @param modelList A list of structural equations created using \code{psem}.
#' @param basis.set An optional list of independence claims.
#' @param direction A \code{vector} of claims defining the specific
#' directionality of independence claims; for use in special cases (see
#' Details).
#' @param interactions whether interactions should be included in independence claims. 
#' Default is FALSE
#' @param conserve Whether the most conservative P-value should be returned;
#' for use in special cases (see Details). Default is FALSE.
#' @param conditioning Whether the conditioning variables should be shown in
#' the summary table. Default is FALSE.
#' @param .progressBar An optional progress bar. Default is TRUE.
#' 
#' @return Returns a \code{data.frame} of independence claims and their
#' significance values.
#' 
#' @author Jon Lefcheck <lefcheckj@@si.edu>, Jarrett Byrnes
#' 
#' @seealso \code{\link{basisSet}}
#' 
#' @references Shipley, Bill. "A new inferential test for path models based on
#' directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.
#' 
#' @export 
#' 
dSep <- function(modelList, basis.set = NULL, direction = NULL, interactions = FALSE, 
                 conserve = FALSE, conditioning = FALSE, .progressBar = TRUE) {
  
  if(is.null(basis.set)) b <- basisSet(modelList, direction, interactions) else b <- basis.set
  
  if(any(duplicated(names(b))) & conserve == FALSE & is.null(direction)) dupOutput(b)
  
  if(length(b) == 0) {
  
    data.frame()
    
  } else {
    
    if(inherits(modelList, "psem")) data <- modelList$data else data <- GetData(modelList)
    
    modelList <- removeData(modelList, formulas = 1)
    
    if(.progressBar == T & length(b) > 0)  pb <- txtProgressBar(min = 0, max = length(b), style = 3)
    
    ret <- do.call(rbind, lapply(1:length(b), function(i)
      
      testBasisSetElements(i, b, modelList, data, conditioning, .progressBar, pb) 
      
      ) )
    
    if(.progressBar == TRUE) close(pb)
    
    rownames(ret) <- NULL
    
    if(conserve == TRUE) {
      
      ret <- do.call(rbind, lapply(unique(names(b)), function(i) {
        
        r <- ret[which(names(b) == i), ]
        
        r[which.min(r[, "P.Value"]), ]
        
      } ) )
      
    }
    
    return(ret)
    
  }
  
}

#' Identify duplicate output
#' 
#' @keywords internal
#' 
dupOutput <- function(b, conserve = FALSE) {
  
  dup <- names(b)[which(duplicated(names(b)))]
  
  if(conserve == FALSE) {
    
    s <- paste("\nNon-linearities detected in the basis set where P-values are not symmetrical.",
               "\nThis can bias the outcome of the tests of directed separation.\n",
               
               "\nOffending independence claims:",
               
               lapply(dup, function(i) {
                 
                 d <- b[names(b) %in% dup]
                 
                 paste(
                   "\n", paste(d[[1]][2], "<-", d[[1]][1]), "*OR*",
                   paste(d[[1]][2], "->", d[[1]][1]), "\n"
                 )
                 
               } ),
               
               "\nOption 1: Specify directionality using argument 'direction = c()' in 'summary'.\n",
               
               "\nOption 2: Remove path from the basis set by specifying as a correlated error using '%~~%' in 'psem'.\n",
               
               "\nOption 3 (recommended): Use argument 'conserve = TRUE' in 'summary' to compute both tests, and return the most conservative P-value.\n"
               
    )
    
    stop(s, call. = FALSE)
    
  }
  
}

#' Evaluate conditional independence claim from the basis set
#' 
#' @keywords internal
#' 
testBasisSetElements <- function(i, b, modelList, data, conditioning, .progressBar, pb) {
  
  formulaList <- lapply(listFormula(modelList, formulas = 1), all.vars_trans, smoothed = TRUE)
  
  bMod <- modelList[[which(sapply(formulaList, function(x) x[1] == b[[i]][2]))]]
  
  # if variable is smoothed and appears in linear model
  
  if(!"gam" %in% class(bMod) & any(grepl("s\\(.*\\)", b[[i]]))) {
    
    bnew <- b[[i]][-2]
    
    bnew <- gsub("(.*)\\,.*", "\\1", gsub("s\\((.*)\\).*", "\\1", bnew)) 
    
    warning("Basis set includes smoothed terms in independence claim: claim is conducted with linear term!", call. = FALSE)
    
  } else bnew <- b[[i]][-2]
  
  if(any(class(bMod) %in% c("lmerMod", "merModLmerTest", "lmerModLmerTest", "glmerMod", "glmmTMB"))) {
    
    bNewMod <- suppressWarnings(
      update(bMod,
             formula(paste(". ~ ", paste(rev(bnew), collapse = " + "), " + ", onlyBars(formula(bMod)))),
             data = data)
    )
    
  } else {
    
    bNewMod <- suppressWarnings(
      update(bMod,
             formula(paste(". ~ ", paste(rev(bnew), collapse = " + "))),
             data = data)
    )
    
  }
  
  ct <- unstdCoefs(bNewMod, data)
  
  ct$Test.Type <- ifelse(is.na(ct$Estimate) | grepl("=", ct$Predictor), "anova", "coef")

  if("gam" %in% class(bNewMod)) {
  
    a <- gsub("(s\\(.*),.*\\)", "\\1", b[[i]][1])
    
    if(any(grepl("s\\(", a))) a <- sapply(a, function(x)
      ifelse(grepl("s\\(", x) & !grepl("\\)", x), paste0(x, ")"), x)) 
  
  } else a <- gsub("(.*)\\,.*", "\\1", gsub("s\\((.*)\\).*", "\\1", b[[i]][1])) 
  
  ct <- ct[which(a == ct$Predictor), , drop = FALSE]
  
  rhs <- paste0(b[[i]][-2], collapse = " + ")
  
  if(conditioning == FALSE) rhs <- paste0(b[[i]][1], " + ...")
  
  ret <- data.frame(Independ.Claim = paste(b[[i]][2], "~", rhs), ct[, c("Test.Type", "DF", "Crit.Value", "P.Value")])
  
  ret <- cbind(ret, isSig(ret[, "P.Value"]))

  names(ret)[ncol(ret)] <- ""
  
  if(.progressBar == TRUE) setTxtProgressBar(pb, i)
  
  return(ret)
  
}

#' Identify duplicate output
#' 
#' @keywords internal
#' 
dupOutput <- function(b, conserve = FALSE) {

  dup <- names(b)[which(duplicated(names(b)))]

  if(conserve == FALSE) {

    s <- paste("\nNon-linearities detected in the basis set where P-values are not symmetrical.",
               "\nThis can bias the outcome of the tests of directed separation.\n",

               "\nOffending independence claims:",

               lapply(dup, function(i) {

                 d <- b[names(b) %in% dup]

                 paste(
                   "\n", paste(d[[1]][2], "<-", d[[1]][1]), "*OR*",
                   paste(d[[1]][2], "->", d[[1]][1]), "\n"
                 )

               } ),

               "\nOption 1: Specify directionality using argument 'direction = c()' in 'summary'.\n",

               "\nOption 2: Remove path from the basis set by specifying as a correlated error using '%~~%' in 'psem'.\n",

               "\nOption 3 (recommended): Use argument 'conserve = TRUE' in 'summary' to compute both tests, and return the most conservative P-value.\n"

               )

    stop(s, call. = FALSE)

  }

}

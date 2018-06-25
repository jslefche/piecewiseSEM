#' Summarizing piecewise structural equation models
#'
#' Returns information necessary to interpret piecewise structural equation
#' models, including tests of directed separation, path coefficients,
#' information criterion values, and R-squared values of individual models.
#'
#' The forthcoming argument \code{groups} splits the analysis based on an optional grouping
#' factor, conducts separate d-sep tests, and reports goodness-of-fit and path
#' coefficients for each submodel. The procedure is approximately similar to a
#' multigroup analysis in traditional variance-covariance SEM. Coming in version 2.1.
#'
#' In cases involving non-normally distributed responses in the independence
#' claims that are modeled using generalized linear models, the significance of
#' the independence claim is not reversable (e.g., the P-value of Y ~ X is not
#' the same as X ~ Y). This is due to the transformation of the response via
#' the link function. In extreme cases, this can bias the goodness-of-fit
#' tests. \code{summary.psem} will issue a warning when this case is present
#' and provide guidance for solutions. One solution is to specify the
#' directionality of the relationship using the \code{direction} argument, e.g.
#' \code{direction = c("X <- Y")}. Another is to run both tests (Y ~ X, X ~ Y)
#' and return the most conservative (i.e., lowest) P-value, which can be
#' toggled using the \code{conserve = TRUE} argument.
#'
#' In some cases, additional claims that were excluded from the basis set can
#' be added back in using the argument \code{add.claims}. These could be, for
#' instance, independence claims among exogenous variables. See Details in
#' \code{\link{basisSet}}.
#'
#' Standardized path coefficients are scaled by standard deviations.
#'
#' @param object a list of structural equations
#' @param ... additional arguments to summary
#' @param direction a vector of claims defining the specific directionality of any independence 
#' claim(s)
#' @param conserve whether the most conservative P-value should be returned (See Details) 
#' Default is FALSE
#' @param conditional whether all conditioning variables should be shown in the table
#' Default is FALSE
#' @param add.claims an optional vector of additional independence claims (P-values) 
#' to be added to the basis set
#' @param intercepts whether intercepts should be included in the coefficient  table
#' Default is FALSE
#' @param standardize whether standardized path coefficients should be reported 
#' Default is "scale"
#' @param standardize.type the type of standardized for non-Gaussian responses: 
#' \code{latent.linear} (default), \code{Mendard.OE}
#' @param .progressBar an optional progress bar. Default is TRUE
#' 
#' @return The function \code{summary.psem} returns a list of summary
#' statistics: \item{dTable}{ A summary table of the tests of directed
#' separation, from \code{\link{dSep}}.  } \item{CStat}{ Fisher's C statistic,
#' degrees of freedom, and significance value based on a Chi-square test.  }
#' \item{IC}{ Information criterion (Akaike, Bayesian, corrected Akaike) as
#' well as degrees of freedom and sample size.  } \item{coefficients}{ A
#' summary table of the path coefficients, from \code{link{coefs}}.  }
#' \item{R2}{ (Pseudo)-R2 values, from \code{\link{rsquared}}.  }
#' @author Jon Lefcheck <jlefcheck@@bigelow.org>
#' @seealso The model fitting function \code{\link{psem}}.
#' @references Shipley, Bill. "A new inferential test for path models based on
#' directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.
#'
#' Shipley, Bill. Cause and correlation in biology: a user's guide to path
#' analysis, structural equations and causal inference. Cambridge University
#' Press, 2002.
#'
#' Shipley, Bill. "Confirmatory path analysis in a generalized multilevel
#' context." Ecology 90.2 (2009): 363-368.
#'
#' Shipley, Bill. "The AIC model selection method applied to path analytic
#' models compared using a d-separation test." Ecology 94.3 (2013): 560-564.
#' 
#' @method summary psem
#'
#' @export
#' 
summary.psem <- function(object, ...,
                         direction = NULL, conserve = FALSE, conditional = FALSE,
                         add.claims = NULL,
                         standardize = "scale", standardize.type = "latent.linear",
                         intercepts = FALSE,
                         .progressBar = TRUE) {

  name <- deparse(substitute(object))

  call <- paste(listFormula(object), collapse = "\n  ")

  dTable <- dSep(object, direction, conserve, conditional, .progressBar)

  Cstat <- fisherC(dTable, add.claims, direction, conserve, conditional, .progressBar)

  IC <- infCrit(object, Cstat, add.claims, direction, conserve, conditional, .progressBar)

  coefficients <- coefs(object, standardize, standardize.type, intercepts)

  R2 <- rsquared(object)

  R2[, which(sapply(R2, is.numeric))] <- round(R2[, which(sapply(R2, is.numeric))], 2)

  if(length(dTable) > 0)

    dTable[, which(sapply(dTable, is.numeric))] <- round(dTable[, which(sapply(dTable, is.numeric))], 4)

  l <- list(name = name, call = call, dTable = dTable, Cstat = Cstat, IC = IC, coefficients = coefficients, R2 = R2)

  class(l) <- "summary.psem"

  l

} 

#' Print summary
#' 
#' @param x an object of class summary.psem
#' @param ... further arguments passed to or from other methods
#' 
#' @method print summary.psem
#' 
#' @export
#' 
print.summary.psem <- function(x, ...) {

  cat("\nStructural Equation Model of", as.character(x$name), "\n")

  cat("\nCall:\n ", x$call)

  cat("\n")

  cat("\n    AIC      BIC")
  cat("\n", as.character(sprintf("%.3f", x$IC[1])), " ", as.character(x$IC[3]))

  cat("\n")

  cat("\n---\nTests of directed separation:\n\n", captureTable(x$dTable))

  cat("\nGlobal goodness-of-fit:\n\n  Fisher's C =", as.character(x$Cstat[1]),
      "with P-value =", as.character(x$Cstat[3]),
      "and on", as.character(x$Cstat[2]), "degrees of freedom")

  cat("\n\n---\nCoefficients:\n\n", captureTable(x$coefficients))

  cat("\n  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05")

  cat("\n\nIndividual R-squared:\n\n", captureTable(x$R2[, c(1, 4:ncol(x$R2))]))

  invisible(x)

}

#' Captures output table
#' 
#' @keywords internal
#' 
captureTable <- function(g) {

  g1 <- capture.output(print(g, row.names = FALSE))
  
  if(all(g1 == "data frame with 0 columns and 0 rows")) 
    
    g1 <- "No independence claims present. Tests of directed separation not possible."

  g1 <- paste0(g1, "\n")

  return(g1)

}

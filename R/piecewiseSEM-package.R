#' Piecewise Structural Equation Modeling
#'
#' Fitting and evaluation of piecewise structural equation models, complete
#' with goodness-of-fit tests, estimates of (standardized) path coefficients,
#' and evaluation of individual model fits (e.g., through R-squared values).
#' Compared with traditional variance-covariance based SEM, piecewise SEM
#' allows for fitting of models to different distributions through GLM and/or
#' hierarchical/nested random structures through (G)LMER. Supported model
#' classes include: \code{lm, glm, gls, Sarlm, lme, glmmPQL, lmerMod, merModLmerTest, 
#' glmerMod, glmmTMB, gam}.
#'
#' \tabular{ll}{ Package: \tab piecewiseSEM\cr Type: \tab Package\cr Version:
#' \tab 2.3.01\cr Date: \tab 2023-06-01\cr Depends: \tab R (>= 4.2.0), car,
#' DiagrammeR, emmeans, igraph, lme4, multcomp, MuMIn, MASS, methods, nlme\cr 
#' License: \tab MIT\cr } The primary functions in the package are \code{\link{psem}} 
#' which unites structural equations in a single model, and \code{summary.psem} can
#' be used on an object of class \code{psem} to provide various summary statistics for 
#' evaluation and interpretation.
#'
#' @name piecewiseSEM-package
#' @aliases piecewiseSEM-package piecewiseSEM
#' @docType package
#' @author Jon Lefcheck <jlefcheck@@umces.edu>
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
#' Shipley, Bill, and Jacob C. Douma. "Generalized AIC and chi‐squared statistics
#'  for path models consistent with directed acyclic graphs." Ecology 101.3 (2020): e02960.
#' 
#' Grace, J.B., Johnson, D.A., Lefcheck, J.S., and Byrnes, J.E.
#' "Standardized Coefficients in Regression and Structural Models with Binary Outcomes."
#' Ecosphere 9(6): e02283.
#'
#' Nakagawa, Shinichi, Paul CD Johnson, and Holger Schielzeth. "The coefficient
#' of determination R2 and intra-class correlation coefficient from generalized
#' linear mixed-effects models revisited and expanded."
#' Journal of the Royal Society Interface 14.134 (2017): 20170213.
#' @keywords package
#'
NULL

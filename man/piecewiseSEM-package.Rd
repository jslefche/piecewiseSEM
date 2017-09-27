\name{piecewiseSEM-package}
\alias{piecewiseSEM-package}
\alias{piecewiseSEM}
\docType{package}
\title{
Piecewise Structural Equation Modeling
}
\description{
Fitting and evaluation of piecewise structural equation models, complete with goodness-of-fit tests, estimates of (standardized) path coefficients, and evaluation of individual model fits (e.g., through R-squared values). Compared with traditional variance-covariance based SEM, piecewise SEM allows for fitting of models to different distributions through GLM and/or hierarchical/nested random structures through (G)LMER. Supported model classes include: \code{lm, glm, gls, pgls, sarlm, lme, glmmPQL, lmerMod, merModLmerTest, glmerMod}.
}
\details{
\tabular{ll}{
Package: \tab piecewiseSEM\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2017-09-27\cr
Depends: \tab R (3.4), nlme, pbkrtest\cr
Suggests: \tab MASS, lme4\cr
License: \tab MIT\cr
}
The primary functions in the package are \code{\link{psem}} which unites structural equations in a single model. \code{summary.psem} can be used on an object of class \code{psem} to provide various summary statistics for evaluation and interpretation.
}
\author{
Jon Lefcheck <jslefche@vims.edu>
}
\references{
Shipley, Bill. "A new inferential test for path models based on directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.

Shipley, Bill. Cause and correlation in biology: a user's guide to path analysis, structural equations and causal inference. Cambridge University Press, 2002.

Shipley, Bill. "Confirmatory path analysis in a generalized multilevel context." Ecology 90.2 (2009): 363-368.

Shipley, Bill. "The AIC model selection method applied to path analytic models compared using a d-separation test." Ecology 94.3 (2013): 560-564.
}
\keyword{ package }
\examples{

}

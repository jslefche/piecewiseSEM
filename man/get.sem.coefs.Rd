\name{get.sem.coefs}
\alias{get.sem.coefs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieves path coefficients for structural equation model.
}
\description{
Extracts unstandardized or standardized (in units of standard deviations of the mean) path coefficients for a piecewise structural equation model (SEM).
}
\usage{
get.sem.coefs(modelList, data, standardized = FALSE, 
corr.errors = NULL)
}
\arguments{
  \item{modelList}{
A list of regressions representing the structural equation model. 
}
\item{data}{
A \code{data.frame} used to construct the structured equations.
}
  \item{standardized}{
Returns standardized coefficients.
}
  \item{corr.errors}{
A vector of variables with correlated errors (separated by "~~").
}
}
\details{
Unstandardized coefficients are extracted directly from model output. Standardized coefficients are scaled by mean and variance using the function \code{scale} on the raw data. If predictors are factors, then they are left in their original units. If models are not gaussian, then the response is left unscaled so the same model parameters can be fit as in the original model.

Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. If two exogenous variables are listed, then the  Spearman correlation between the two is given, otherwise the correlation between the partial residuals is returned. P-values for correlated errors are derived from a t-distribution using the transformation in Shipley, 2000 (p.76).
}
\value{
Returns a \code{data.frame} with the causal path, the coefficient estimate, the standard error, and p-value.
}
\references{
Shipley, Bill. Cause and correlation in biology: a user's guide to path analysis, structural equations and causal inference. Cambridge University Press, 2002.
}
\author{
Jon Lefcheck
}
\note{
}
\seealso{
\code{\link{lm.beta}}, \code{\link{scale}}, \code{\link{summary}}
}
\examples{
# Fit model from Shipley (2009)
data(shipley2009) 

# Load model packages
library(lmerTest)
library(nlme)

# Create list of models 
shipley2009.modlist = list(
  lme(DD~lat, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  lme(Date~DD, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  lme(Growth~Date, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  glmer(Live~Growth+(1|site)+(1|tree), 
  family=binomial(link = "logit"), data = shipley2009) )
  
# Get path coefficients
get.sem.coefs(shipley2009.modlist, shipley2009)
get.sem.coefs(shipley2009.modlist, shipley2009, standardized = TRUE)
}
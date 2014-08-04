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
get.sem.coefs(modelList, standardized = FALSE)
}
\arguments{
  \item{modelList}{
A list of regressions representing the structural equation model. 
}
  \item{standardized}{
Returns standardized coefficients.
}
}
\details{
Unstandardized coefficients are extracted directly from model output.

Standardized coefficients are scaled by mean and variance using the function \code{scale} on the raw data. If predictors are factors, then they are left in their original units. If models are not gaussian, then the response is left unscaled so the same model parameters can be fit as in the original model.
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
data(shipley) # Data from Shipley 2009

# Load model packages
library(lmerTest)
library(nlme)

# Create list of models 
Shipley.modlist = list(
  lme(DD~lat, random = ~1|site/tree, na.action = na.omit, 
  data = Shipley),
  lme(Date~DD, random = ~1|site/tree, na.action = na.omit, 
  data = Shipley),
  lme(Growth~Date, random = ~1|site/tree, na.action = na.omit, 
  data = Shipley),
  glmer(Live~Growth+(1|site)+(1|tree), 
  family=binomial(link = "logit"), data = Shipley) )
  
# Get path coefficients
get.sem.coefs(Shipley.modlist)
}
\name{get.fisher.c}
\alias{get.fisher.c}
\title{
Goodness-of-fit for piecewise SEM.
}
\description{
Evaluates independence claims for a piecewise structural equation model (SEM), and calculates Fisher's C and associated p-value to evaluate model fit.
}
\usage{
get.fisher.c(modelList, corr.errors = NULL, add.vars = NULL, adjust.p = FALSE, 
basis.set = NULL, pvalues.df = NULL, .progressBar = TRUE)
}
\arguments{
  \item{modelList}{
A \code{list} of regressions representing the structural equation model. 
}
  \item{corr.errors}{
A vector of variables with correlated errors (separated by "~~").
}
  \item{add.vars}{
A vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
  \item{adjust.p}{
Whether p-values degrees of freedom should be adjusted.
}
  \item{basis.set}{
Provide an optional basis set.
}
  \item{pvalues.df}{
An optional \code{data.frame} corresponding to p-values for independence claims.
}
  \item{.progressBar}{
Enable optional text progress bar.
}
}
\details{
}
\value{
Returns a vector with the first entry corresponding to Fisher's C statistic, and the second corresponding to the associated p-value (derived from a Chi-squared distribution). 
}
\references{
Shipley, Bill. "Confirmatory path analysis in a generalized multilevel context." Ecology 90.2 (2009): 363-368.
}
\author{
Jon Lefcheck
}
\note{
}
\seealso{
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

# Get goodness-of-fit statistics
get.fisher.c(modelList)
}
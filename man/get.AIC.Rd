\name{get.aic}
\alias{get.aic}
\title{
Extracts AIC scores for a piecewise structural equation model.
}
\description{
Extracts the AIC and AICc (corrected for small sample size) values from a piecewise structural equation model (SEM).
}
\usage{
get.aic(modelList, data, corr.errors = NULL, 
add.vars = NULL, grouping.var = NULL, 
top.level.vars = NULL, adjust.p = FALSE,
basis.set = NULL, pvalues.df = NULL, 
.progressBar = TRUE)  
}
\arguments{
  \item{modelList}{
A \code{list} of regressions representing the structural equation model. 
}
\item{data}{
A \code{data.frame} used to construct the structured equations.
}
  \item{corr.errors}{
A vector of variables with correlated errors (separated by "~~").
}
  \item{add.vars}{
A vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
\item{grouping.var}{
An optional variable that represents the highest level of a hierarchical dataset.
}
\item{top.level.vars}{
An optional \code{vector} of variables that are identical at the highest level of the hierarchy.
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
This function calculates AIC and AICc (corrected for small sample sizes) values for a piecewise structural equation model (SEM). For more information, see Shipley 2013.
}
\value{
Returns a vector where the first entry is the AIC score, and the second is the AICc score, and the third is the likelihood degrees of freedom (K).
}
\references{
Shipley, Bill. "The AIC model selection method applied to path analytic models compared using a d-separation test." Ecology 94.3 (2013): 560-564.
}
\author{
Jon Lefcheck
}
\note{
}
\seealso{
}
\examples{
data(Shipley) # Data from Shipley 2009

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
  
# Get AIC scores
get.aic(Shipley.modlist, Shipley)
}
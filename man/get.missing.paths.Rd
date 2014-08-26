\name{get.missing.paths}
\alias{get.missing.paths}
\title{
Evaluate independence claims for piecewise SEM.
}
\description{
Identifies missing paths from a piecewise SEM, fits models, extracts path p-values and returns in a \code{data.frame}.
}
\usage{
get.missing.paths(modelList, adjust.p = FALSE, .progressBar = FALSE, 
basis.set = NULL, add.vars = NULL, corr.errors = NULL)
}
\arguments{
  \item{modelList}{
A \code{list} of regressions representing the structural equation model. 
}
  \item{adjust.p}{
Whether p-values degrees of freedom should be adjusted.
}
  \item{.progressBar}{
Enable optional text progress bar.
}
  \item{basis.set}{
Provide an optional basis set.
}
  \item{add.vars}{
A vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
  \item{corr.errors}{
A vector of variables with correlated errors (separated by "~~").
}
}
\details{
This function takes a model list (and optional basis set) and evaluates all conditional independence claims by constructing regressions, returning the claims, the variables upon which they are conditional, and associated p-values in a \code{data.frame}.

Models are constructed based on how the variables are treated as in the model list. For example, if the indepedence claim includes a binary variable that is fit to a binomial distribution using an identity link, the function will evaluate the claim using the same parameters. Similarly, for linear mixed effects models construted in lme4 or nlme, varying slopes and intercepts are treated as in the model list. For example, if a variable is modeled with both a random slope and intercept in any model in the model list, that variable will be modeled with a random slope and intercept when evaluating all independence claims in which it appears. If slopes and intercepts vary for multiple variables, they will appear as such, even if they are conditional.

For linear mixed effects models, p-values can be adjusted to accommodate the full model degrees of freedom using the argument \code{p.adjust = TRUE}. For more information, see Shipley 2013.

Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{get.sem.coefs}.
}
\value{
Returns a \code{data.frame} with the first column representing the independence claim as a directed arrow between the two variables, the second column representing a vector of variables upon which the claim is conditional, and the third column the p-value from a regression modeling the missing path.
}
\references{
Shipley, Bill. "Confirmatory path analysis in a generalized multilevel context." Ecology 90.2 (2009): 363-368.

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

# Evaluate independence claims
get.missing.paths(modelList, .progressBar = T)
}
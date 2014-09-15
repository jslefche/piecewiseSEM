\name{get.sem.fit}
\alias{get.sem.fit}
\title{
Shipley's goodness-of-fit tests for piecewise SEM.
}
\description{
Tests independence claims and calculates Fisher's C statistic and associated p-value, and AIC and AICc, for piecewise structural equation model (SEM).
}
\usage{
get.sem.fit = function(modelList, data, 
corr.errors = NULL, add.vars = NULL, adjust.p = FALSE, 
basis.set = NULL, pvalues.df = NULL, .progressBar = TRUE)
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
Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{get.sem.coefs}.

The argument \code{add.vars} requires a vector of character strings corresponding to column names in the dataset used to construct the models in \code{modelList}. This is useful if comparing nested SEMs where one wishes to account for additional variables whose independence claims should be evaluated, but which do not have any hypothesized paths in the current SEM. The default is \code{NULL} (there is no additional independence claims that do not appear in the model list).
For linear mixed effects models, p-values can be adjusted to accommodate the full model degrees of freedom using the argument \code{p.adjust = TRUE}. For more information, see Shipley 2013.

Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{get.sem.coefs}.
}
\value{
Returns a \code{list} with the following:
\item{missing.paths}{A \code{data.frame} where the first column are the independence claims, the second column is a vector of variables upon which the claims conditional, and the third column are the p-values for a given claim. }
\item{fishers.c}{ A vector where the first item is Fisher's C statistic, and the second is the associated p-value (derived from a Chi-squared distribution) }
\item{AIC}{ A vector where the first item is the AIC score of the model, and the second is the AICc score (corrected for small sample size). }
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
\code{\link{dag.updated}}, \code{\link{get.missing.paths}}, \code{\link{get.fisher.c}}, \code{\link{get.aic}}
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

# Conduct goodness-of-fit tests
get.sem.fit(Shipley.modlist, Shipley)
}
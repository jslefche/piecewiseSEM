\name{get.sem.fit.multilevel}
\alias{get.sem.fit.multilevel}
\title{
Shipley's goodness-of-fit tests for piecewise SEM.
}
\description{
Tests independence claims and calculates Fisher's C statistic and associated p-value, and AIC and AICc, for a multilevel piecewise structural equation model (SEM).
}
\usage{
get.sem.fit.multilevel = function(modelList, data, 
grouping.var = NULL, top.level.vars = NULL,
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
\item{grouping.var}{
An optional variable that represents the highest level of a hierarchical dataset.
}
\item{top.level.vars}{
An optional \code{vector} of variables that are identical at the highest level of the hierarchy.
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
Replicates the functionality of \code{get.sem.fit} but allows the user to optionally specify a top level variables for a hierarchical dataset, which will be summarized based on the grouping variable. For example, consider a two-level hierarchy, where variables at the top level have identical values for each level of the grouping variable. If the response is a top level variable (is identically replicated for the grouping variable), this function takes a mean of the lower level variables for each level of the grouping variable, then runs the test of d-separation. If the response is fully replicated (occurs at the lower level), then no aggregation occurs.
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

}
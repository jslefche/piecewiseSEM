\name{get.sem.fit}
\alias{get.sem.fit}
\title{
Shipley's goodness-of-fit tests for piecewise SEM.
}
\description{
Tests independence claims and calculates Fisher's C statistic and associated p-value, and AIC and AICc, for a piecewise structural equation model (SEM).
}
\usage{
get.sem.fit = function(modelList, data, 
corr.errors = NULL, add.vars = NULL, adjust.p = FALSE, 
basis.set = NULL, pvalues.df = NULL, 
disp.conditional = FALSE, model.control = NULL,
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
  \item{disp.conditional}{
Display conditional variables in reporting d-sep test results.
}
  \item{model.control}{
A \code{list} of model control arguments to be passed to d-sep models.
}
  \item{.progressBar}{
Enable optional text progress bar.
}
}
\details{
Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{get.sem.coefs}.

The argument \code{add.vars} requires a vector of character strings corresponding to column names in the dataset used to construct the models in \code{modelList}. This is useful if comparing nested SEMs where one wishes to account for additional variables whose independence claims should be evaluated, but which do not have any hypothesized paths in the current SEM. The default is \code{NULL} (there is no additional independence claims that do not appear in the model list).
For linear mixed effects models, p-values can be adjusted to accommodate the full model degrees of freedom using the argument \code{p.adjust = TRUE}. For more information, see Shipley 2013.

\code{top.level.vars} allows the user to optionally specify one or more top level variables for a hierarchical dataset, which will be summarized based on one or more grouping variables specified by \code{grouping.vars}. For example, consider a two-level hierarchy, where variables at the top level have identical values for each level of the grouping variable. If the response is a top level variable (is identically replicated for the grouping variable), this function takes a mean of the lower level variables for each level of the grouping variable, then runs the test of d-separation. If the response is fully replicated (occurs at the lower level), then no aggregation occurs.

Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{get.sem.coefs}.

The model controls in \code{model.control} will be passed to any subsequent d-sep models that fall under that model class. So, if a response is fit using \code{lme} and the \code{list} provided in \code{model.control} contains \code{lmeControl}, then all d-sep models including that response will inherit those arguments from \code{lmeControl}.
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

# Conduct goodness-of-fit tests
get.sem.fitshipley2009, shipley2009)

# Add new variable with no hypothesized links and test for independence
set.seed(1)
add.var = rnorm(nrow(shipley2009), 50, 20)

get.sem.fit(shipley2009.modlist, shipley2009, add.vars = c("add.var"))

# Add new variable grouped at site level
set.seed(2)
shipley2009$site.var = rep(rnorm(20, 10, 5), each = 95)

shipley2009.modlist2 = append(shipley2009, list(
  lm(site.var~Growth, data = aggregate(shipley2009, by = list(site = Shipley[,"site"]), mean, na.rm = T))) )

get.sem.fit(shipley2009.modlist2, shipley2009, grouping.vars = c("site"), top.level.vars = c("site.var"))

# Fit model from Shipley (2013)
data(shipley2013)

shipley2013.modlist = list(
  lme(x2~x1, random = ~x1 | species, data = shipley2013),
  lme(x3~x2, random = ~x2 | species, data = shipley2013),
  lme(x4~x2, random = ~x2 | species, data = shipley2013),
  lme(x5~x3+x4, random = ~x3+x4 | species, data = shipley2013) )

## NOT RUN ##
# get.sem.fit(shipley2013.modlist, shipley2013) 
# Convergence error!
## NOT RUN ##

# Specify model controls by switching to old optimizer
get.sem.fit(shipley2013.modlist, shipley2013, model.control = list(lmeControl(opt = "optim")))

}
\name{get.missing.paths}
\alias{get.missing.paths}
\title{
Evaluate independence claims for piecewise SEM
}
\description{
Identifies missing paths from a piecewise SEM, fits models, extracts path p-values and returns in a \code{data.frame}.
}
\usage{
get.missing.paths(modelList, data, corr.errors, add.vars, grouping.var, 
                  top.level.vars, adjust.p = FALSE, basis.set, 
                  disp.conditional = FALSE, model.control, sig = 3, 
                  .progressBar = TRUE)
}
\arguments{
  \item{modelList}{
a \code{list} of regressions representing the structural equation model. 
}
\item{data}{
a \code{data.frame} used to construct the structured equations.
}
  \item{corr.errors}{
a vector of variables with correlated errors (separated by "~~").
}
  \item{add.vars}{
a vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
\item{grouping.var}{
an optional variable that represents the highest level of a hierarchical dataset.
}
\item{top.level.vars}{
an optional \code{vector} of variables that are identical at the highest level of the hierarchy.
}
  \item{adjust.p}{
whether p-values degrees of freedom should be adjusted (see below). Default is \code{FALSE}.
}
  \item{basis.set}{
provide an optional basis set.
}
  \item{disp.conditional}{
display conditional variables in reporting d-sep test results. Default is \code{FALSE}.
}
  \item{model.control}{
a \code{list} of model control arguments to be passed to d-sep models.
}
\item{sig}{
number of significant digits to round output. Default is \code{3}.
}
  \item{.progressBar}{
enable optional text progress bar. Default is \code{TRUE}.
}
}
\details{
This function takes a model list (and optional basis set) and evaluates all conditional independence claims by constructing regressions, returning the claims, the variables upon which they are conditional, and associated p-values in a \code{data.frame}.

Models are constructed based on how the variables are treated as in the model list. For example, if the indepedence claim includes a binary variable that is fit to a binomial distribution using an identity link, the function will evaluate the any claims using the same parameters. Similarly, for linear mixed effects models construted in lme4 or nlme, varying slopes and intercepts are treated as in the model list. For example, if a variable is modeled with both a random slope and intercept in any model in the model list, that variable will be modeled with a random slope and intercept when evaluating all independence claims in which it appears. If slopes and intercepts vary for multiple variables, they will appear as such, even if they are conditional.

Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{get.sem.coefs}.

The argument \code{add.vars} requires a vector of character strings corresponding to column names in the dataset used to construct the models in \code{modelList}. This is useful if comparing nested SEMs where one wishes to account for additional variables whose independence claims should be evaluated, but which do not have any hypothesized paths in the current SEM. The default is \code{NULL} (there is no additional independence claims that do not appear in the model list).

For linear mixed effects models, p-values can be adjusted to accommodate the full model degrees of freedom using the argument \code{p.adjust = TRUE}. For more information, see Shipley 2013.
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

# Evaluate independence claims
get.missing.paths(shipley2009.modlist, shipley2009)
}
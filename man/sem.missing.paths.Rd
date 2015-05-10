\name{sem.missing.paths}

\alias{sem.missing.paths}

\title{
Evaluate independence claims for piecewise SEM
}

\description{
Identifies missing paths from a piecewise SEM, fits models, extracts path p-values and returns in a \code{data.frame}.
}

\usage{
sem.missing.paths(modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, 
  filter.exog, adjust.p, basis.set, model.control, .progressBar) 
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

\item{grouping.vars}{
an optional variable that represents the highest level of a hierarchical dataset.
}

\item{top.level.vars}{
an optional \code{vector} of variables that are identical at the highest level of the hierarchy.
}

\item{filter.exog}{
filters relationships between exogenous variables from the basis set. Default is \code{TRUE}.
}

\item{adjust.p}{
whether p-values degrees of freedom should be adjusted (see below). Default is \code{FALSE}.
}

\item{basis.set}{
provide an optional basis set.
}

\item{model.control}{
a \code{list} of model control arguments to be passed to d-sep models.
}

\item{.progressBar}{
enable optional text progress bar. Default is \code{TRUE}.
}

}

\details{
This function takes a model list (and optional basis set) and evaluates all conditional independence claims by constructing regressions, returning the claims, the variables upon which they are conditional, and associated p-values in a \code{data.frame}.

Models are constructed based on how the variables are treated as in the model list. For example, if the indepedence claim includes a binary variable that is fit to a binomial distribution using an identity link, the function will evaluate the any claims using the same parameters. Similarly, for linear mixed effects models construted in lme4 or nlme, varying slopes and intercepts are treated as in the model list. For example, if a variable is modeled with both a random slope and intercept in any model in the model list, that variable will be modeled with a random slope and intercept when evaluating all independence claims in which it appears. If slopes and intercepts vary for multiple variables, they will appear as such, even if they are conditional.
}

\value{
Returns a \code{data.frame} where the first column is the independence claim, and the second through sixth columns the model estimates corresponding to the response variable in the independence claim.
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

# Load model packages
library(lmerTest)
library(nlme)

# Load data from Shipley (2009)
data(shipley2009) 

# Create list of models 
shipley2009.modlist = list(

  lme(DD~lat, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  
  lme(Date~DD, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  
  lme(Growth~Date, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  
  glmer(Live~Growth+(1|site)+(1|tree), 
  family=binomial(link = "logit"), data = shipley2009) 
  
  )
  
# Evaluate independence claims
sem.missing.paths(shipley2009.modlist, shipley2009)

}
\name{sem.aic}

\alias{sem.aic}

\title{
Extracts AIC scores for piecewise SEM
}

\description{
Extracts the AIC and AICc (corrected for small sample size) values from a piecewise structural equation model (SEM).
}

\usage{
sem.aic(modelList, data, corr.errors, add.vars, grouping.vars, top.level.vars, 
  filter.exog, adjust.p, basis.set, pvalues.df, model.control, .progressBar) 
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

\item{pvalues.df}{
an optional \code{data.frame} corresponding to p-values for independence claims.
}

\item{model.control}{
a \code{list} of model control arguments to be passed to d-sep models.
}

\item{.progressBar}{
enable optional text progress bar. Default is \code{TRUE}.
}
}

\details{
This function calculates AIC and AICc (corrected for small sample sizes) values for a piecewise structural equation model (SEM).

For linear mixed effects models, p-values can be adjusted to accommodate the full model degrees of freedom using the argument \code{p.adjust = TRUE}. For more information, see Shipley 2013.
}

\value{
Returns a \code{data.frame} where the first entry is the AIC score, and the second is the AICc score, and the third is the likelihood degrees of freedom (K).
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
  
# Get AIC and AICc values for the SEM
sem.aic(shipley2009.modlist, shipley2009)

}
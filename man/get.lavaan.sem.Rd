\name{get.lavaan.sem}
\alias{get.lavaan.sem}
\title{
Constructs variance-covariance based SEM
}
\description{
Estimates variance-covariance based (traditional) structural equation model (SEM) using the \code{lavaan} package.
}
\usage{
get.lavaan.sem(modelList, data, corr.errors, add.vars, ...)
}
\arguments{
  \item{modelList}{
a list of regressions representing the structural equation model. 
}
  \item{data}{
a \code{data.frame} containing the full dataset for the SEM.
}
  \item{corr.errors}{
a vector of variables with correlated errors (separated by "~~").
}
  \item{add.vars}{
a vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
  \item{...}{
additional arguments passed to \code{lavaan}.
}
}
\details{
}
\value{
Returns a model object of class \code{lavaan}.
}
\references{
Grace, James B. Structural equation modeling and natural systems. Cambridge University Press, 2006.
}
\author{
Jon Lefcheck
}
\note{
}

\seealso{
\code{\link[lavaan]{sem}}
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
  
# Get variance-covariance based SEM
summary(get.lavaan.sem(shipley2009.modlist, shipley2009))
}
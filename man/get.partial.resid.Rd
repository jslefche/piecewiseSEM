\name{get.partial.resid}
\alias{get.partial.resid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates partial residuals for two variables.
}
\description{
Extracts partial residuals for \code{y ~ x | Z}, where \code{Z} represents all other variables present in a structured equation upon which \code{y} and \code{x} are conditional.
}
\usage{
get.partial.resid(y ~ x, modelList, model.control = NULL)
}
\arguments{
   \item{.formula}{
The partial residuals of \code{y ~ x}.
}
 \item{modelList}{
A list of regressions representing the structural equation model. 
}
  \item{model.control}{
A \code{list} of model control arguments to be passed to partial models.
}
  \item{plot}{
Whether the partial plot of \code{y ~ x | Z} should be returned. Default is \code{TRUE}.
}
  \item{regr}{
Whether the partial regression of \code{y ~ x | Z} should also be plotted. Default is \code{TRUE}.
}
}
\details{

}
\value{
Returns a \code{data.frame} where the first column are the partial residuals of \code{y ~ x | Z} , and the second column is the partial residuals of \code{x ~ Z | y}. If \code{plot == TRUE} then a plot of \code{y ~ x | Z}. If \code{regr == TRUE} then the plot includes the linear regression of \code{y ~ x | Z} .
}
\references{
Shipley, Bill. Cause and correlation in biology: a user's guide to path analysis, structural equations and causal inference. Cambridge University Press, 2002.
}
\author{
Jon Lefcheck
}
\note{
}
\seealso{

}
\examples{
# Fit model from Shipley (2013)
data(shipley2013)

# Load model packages
library(nlme)

# Create list of models
shipley2013.modlist = list(
  lme(x2~x1, random = ~x1 | species, data = shipley2013),
  lme(x3~x2, random = ~x2 | species, data = shipley2013),
  lme(x4~x2, random = ~x2 | species, data = shipley2013),
  lme(x5~x3+x4, random = ~x3+x4 | species, data = shipley2013) )

# Get partial residuals of x3 on x5 conditional on x4
resids.df = get.partial.resid(x5 ~ x3, shipley2013.modlist, list(lmeControl(opt = "optim")))

# Also returns raw residuals values for plotting in other packages
head(resids.df)
}
\name{get.partial.resid}
\alias{get.partial.resid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates partial residuals for two variables.
}
\description{
Extracts partial residuals for y ~ x | Z, where Z represents all other variables present in a structured equation upon which y and x are conditional.
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
}
\details{

}
\value{
Returns a \code{data.frame} where the first column are the partial residuals of y ~ Z | x, and the second column is the partial residuals of x ~ Z | y.
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
resids.df = get.partial.resid(x5 ~ x3, shipley2013.modlist, shipley2013, list(lmeControl(opt = "optim")))

# Plot partial residuals against one another to look for partial dependency
plot(resids.df[ ,1] ~ resids.df [,2])
}
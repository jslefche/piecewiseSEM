\name{get.model.fit}
\alias{get.model.fit}
\title{
Goodness-of-fit statistics for component models.
}
\description{
Returns (pseudo)-R^2 and AIC values for component models in structural equation model (SEM).
}
\usage{
get.sem.fit = function(modelList)
}
\arguments{
   \item{modelList}{
A \code{list} of regressions representing the structural equation model. 
}

}
\details{
Returns goodness-of-fit statistics for generalized linear (mixed) models, including (marginal and condition) R^2 and Akaike Information Criterion (AIC) values. 
}
\value{
Returns a \code{data.frame} with the model class, the family, the link function,  the marginal R^2 (based on fixed effects only) and/or conditional R^2 (based on fixed and random effects, if present), and the AIC score.
}
\references{
Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method for obtaining R2 from 
generalized linear mixed‐effects models." Methods in Ecology and Evolution 4.2 (2013): 133-142.

Johnson, Paul C.D. "Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models." Methods in
Ecology and Evolution.
}
\author{
Jon Lefcheck & Juan Casallas
}
\note{
This function sources the \code{rsquared.glmm} function found here: https://github.com/jslefche/rsquared.glmm, and thus may be periodically updated independently of this package.
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

# Return model fit statistics
get.model.fits(shipley2009.modlist)
}
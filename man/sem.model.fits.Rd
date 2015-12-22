\name{sem.model.fits}

\alias{sem.model.fits}

\title{
Goodness-of-fit statistics for linear models
}

\description{
Returns (pseudo)-R^2 and AIC values for component models in structural equation model (SEM).
}

\usage{
sem.model.fits(modelList, aicc = FALSE)
}

\arguments{
\item{modelList}{
a \code{list} of regressions representing the structural equation model. 
}

\item{aicc}{
whether AIC corrected for small sample size (AICc) should be returned. Default is \code{FALSE}. 
}

}
\details{
Returns goodness-of-fit statistics for generalized linear (mixed) models, including (marginal and condition) R^2 and Akaike Information Criterion (AIC(c)) values. 
}

\value{
Returns a \code{data.frame} with the model class, the family, the link function,  the marginal R^2 (based on fixed effects only) and/or conditional R^2 (based on fixed and random effects, if present), and the AIC(c) score (based on ML).
}

\references{
Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method for obtaining R2 from 
generalized linear mixed-effects models." Methods in Ecology and Evolution 4.2 (2013): 133-142.

Johnson, Paul C.D. "Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models." Methods in
Ecology and Evolution.
}

\author{
Jon Lefcheck & Juan Casallas
}

\note{
This function sources the \code{rsquared.glmm} function found here: https://github.com/jslefche/rsquared.glmm, and thus may be periodically updated independently of this package.
}

\examples{
# Load data from Shipley (2009)
data(shipley2009) 

# Get R2 for linear model
lm.mod = lm(DD ~ lat, data = shipley2009)
sem.model.fits(lm.mod)

# Get R2 for generalized linear model
glm.mod = glm(Live ~ Growth, family = "binomial", data = shipley2009)
sem.model.fits(glm.mod)

# Get R2 for generalized least-squares model
library(nlme) 

gls.mod = gls(DD ~ lat, na.action = na.omit, data = shipley2009)
sem.model.fits(gls.mod)

# Can supply the models as a list
# Use lm and gls -- should produce very similar R2s
sem.model.fits(list(lm.mod, gls.mod))

# Get R2 for linear mixed effects model (nlme)
lme.mod = lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit, data = shipley2009)
sem.model.fits(lme.mod)

# Get R2 for linear mixed effects model (lme4)
library(lme4)

lmer.mod = lmer(DD ~ lat + (1|site/tree), data = shipley2009)
sem.model.fits(lmer.mod)

# Get R2 for generalized linear mixed effects model (lme4)
glmer.mod = glmer(Live ~ Growth + (1|site/tree), family = "binomial", data = shipley2009)
sem.model.fits(glmer.mod)

### SEM EXAMPLE ###

# Create list of models 
shipley2009.modlist = list(
  
  lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Date ~ DD, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Growth ~ Date, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  glmer(Live ~ Growth+(1|site)+(1|tree), 
        family=binomial(link = "logit"), data = shipley2009) 
  
)

# Return model fit statistics
sem.model.fits(shipley2009.modlist)
}
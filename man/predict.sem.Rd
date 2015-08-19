\name{predict.sem}

\alias{predict.sem}

\title{
Returns model predictions for SEM
}

\description{
Returns predictions for responses in the list of structured equations.
}

\usage{
predict.sem(modelList, newdata, ...)
}

\arguments{

\item{modelList}{
a single model or list of regressions representing the structural equation model. 
}

\item{newdata}{
a \code{data.frame} of predictors used to generate the model predictions.
}

\item{sefit}{
whether standard errors of predictions should be returned. Default is \code{FALSE}.
}

\item{...}{
additional arguments passed to \code{predict}.
}

}

\details{
Mixed model predictions includes only fixed effects and *not* random effects (default is \code{level = 0} or \code{re.form = 0}). This can be changed by passing additional arguments from \code{predict.lme} and \code{predict.merMod}.

If \code{sefit = TRUE} for mixed models, then standard errors on predictions are estimated using fixed effects *only*. See explanation here: http://glmm.wikidot.com/faq.
}

\value{
Returns a \code{data.frame} containing the new data and the predicted responses.
}

\references{
}

\author{
Jon Lefcheck
}

\note{
}

\seealso{
\code{\link{predict}}, \code{\link{predict.lme}}, \code{\link{predict.merMod}}
}

\examples{
# Load example data
data(shipley2009)

# Load model packages
library(lmerTest)
library(nlme)

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

# Create new data for predictions
shipley2009.new = data.frame(
  
  lat = rnorm(length(shipley2009$lat), 
    mean(shipley2009$lat, na.rm = TRUE), 
    sd(shipley2009$lat, na.rm = TRUE)),
  
  DD = rnorm(length(shipley2009$DD), 
    mean(shipley2009$DD, na.rm = TRUE), 
    sd(shipley2009$DD, na.rm = TRUE)),
    
  Date = rnorm(length(shipley2009$Date), 
    mean(shipley2009$Date, na.rm = TRUE), 
    sd(shipley2009$Date, na.rm = TRUE)),
    
  Growth = rnorm(length(shipley2009$Growth), 
    mean(shipley2009$Growth, na.rm = TRUE), 
    sd(shipley2009$Growth, na.rm = TRUE))
  
)

# Generate predictions
head(predict.sem(shipley2009.modlist, shipley2009.new))

# Generate predictions with standard errors (based on fixed effects only)
head(predict.sem(shipley2009.modlist, shipley2009.new, sefit = TRUE))

###

# Load data from Shipley (2013)
data(shipley2013) 

shipley2013.modlist = list(

  lme(x2~x1, random = ~x1 | species, data = shipley2013),
  
  lme(x3~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x4~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x5~x3+x4, random = ~x3+x4 | species, data = shipley2013)
  
  )

# Generate new data
shipley2013.new = data.frame(

  x1 = rnorm(length(shipley2013$x1), 
    mean(shipley2013$x1), 
    sd(shipley2013$x1)),
  
  x2 = rnorm(length(shipley2013$x2), 
    mean(shipley2013$x2), 
    sd(shipley2013$x2)),
  
  x3 = rnorm(length(shipley2013$x3), 
    mean(shipley2013$x3), 
    sd(shipley2013$x3)),
  
  x4 = rnorm(length(shipley2013$x4), 
    mean(shipley2013$x4), 
    sd(shipley2013$x4))
  
)

# Generate predictions 
head(predict.sem(shipley2013.modlist, shipley2013.new))
}
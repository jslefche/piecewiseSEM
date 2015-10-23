\name{sem.predict}

\alias{sem.predict}

\title{
Returns model predictions for SEM
}

\description{
Returns predictions for responses in the list of structured equations.
}

\usage{
sem.predict(object, newdata, sefit = FALSE, ...)
}

\arguments{

\item{object}{
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
Returns a \code{data.frame} containing the new data and the predicted responses based on fixed effects only for the innermost level of grouping (i.e., population residuals) corresponding to \code{re.form = NA} for \code{lme4} models and \code{level = Q} for \code{nlme} models.
}

\author{
Jon Lefcheck
}

\seealso{
\code{\link{predict}}, \code{\link[nlme]{predict.lme}}, \code{\link[lme4]{predict.merMod}}
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
  DD = seq(min(shipley2009$DD, na.rm = TRUE), 
           max(shipley2009$DD, na.rm = TRUE), 
           by = 0.01)
)

# Generate predictions
shipley2009.new.pred = sem.predict(shipley2009.modlist, shipley2009.new)
head(shipley2009.new.pred)

# Plot predicted fit
plot(shipley2009$Date ~ shipley2009$DD, col = "grey60")
lines(shipley2009.new.pred$Date.fit ~ shipley2009.new.pred$DD, lwd = 2, col = "red")

# Generate predictions with standard errors (based on fixed effects only)
shipley2009.new.pred = sem.predict(shipley2009.modlist, shipley2009.new, sefit = TRUE)

# Add 95% confidence intervals
lines(shipley2009.new.pred$DD, 
      shipley2009.new.pred$Date.fit + 2 * shipley2009.new.pred$Date.se.fit, 
      lwd = 1.5, lty = 2, col = "red")

lines(shipley2009.new.pred$DD, 
      shipley2009.new.pred$Date.fit - 2 * shipley2009.new.pred$Date.se.fit,
      lwd = 1.5, lty = 2, col = "red")
}
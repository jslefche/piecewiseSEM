\name{partial.resid}

\alias{partial.resid}

\title{
Calculates partial residuals for two variables
}

\description{
Extracts partial residuals for \code{y ~ x | Z}, where \code{Z} represents all other variables present in a structured equation upon which \code{x} is conditional.
}

\usage{
partial.resid(.formula = y ~ x, modelList, data, model.control, return.data.frame = FALSE,
plotit = TRUE, plotreg = TRUE, plotCI = TRUE)
}

\arguments{

\item{.formula}{
a \code{formula} specifying the partial residuals to investigate.
}

\item{modelList}{
a single model or list of regressions representing the structural equation model. 
}

\item{data}{
a \code{data.frame} used to construct the model.
}

\item{model.control}{
a \code{list} of model control arguments to be passed to the partial residual models.
}

\item{return.data.frame}{
whether a \code{data.frame} of the partial residuals should be returned. Default is \code{TRUE}.
}

\item{plotit}{
whether the partial plot of \code{y ~ x | Z} should be returned. Default is \code{TRUE}.
}

\item{plotreg}{
whether the partial regression of \code{y ~ x | Z} should also be plotted. Default is \code{TRUE}.
}

\item{plotCI}{
whether the confidence intervals of the partial regression of \code{y ~ x | Z} should also be plotted. Default is \code{TRUE}.
}

}

\details{

}
\value{
Returns a \code{data.frame} where the first column are the partial residuals of \code{y ~ x | Z}, and the second column is the partial residuals of \code{x ~ Z | y}. 

If \code{plotit = TRUE}, then the function also returns a plot of \code{y ~ x | Z}. 

If \code{plotreg = TRUE}, then the plot includes the linear regression of \code{y ~ x | Z}.
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
# Create example data
set.seed(1)

example.data = data.frame(
  y = rnorm(100, 0, 1),
  x1 = rnorm(100, 10, 50),
  random = letters[1:5]
)

example.data$x2 = example.data$y + runif(100, 0, 5)

example.data$x3 = example.data$x2 + runif(100, 0, 5)

# Run regular linear model using lm()
lm.model = lm(y ~ x1 + x2 + x3, example.data)

partial.resid(y ~ x3, lm.model, example.data)

# Works with interactions too
lm.model2 = lm(y ~ x1 * x2 + x3, example.data)

partial.resid(y ~ x1 * x2, lm.model2, example.data, return.data.frame = FALSE)

# Run generalized least squared regression
library(nlme)

gls.model = gls(y ~ x1 + x2 + x3, example.data)

partial.resid(y ~ x3, gls.model, example.data, return.data.frame = FALSE)

# Run mixed effects model using lme()
lme.model = lme(y ~ x1 + x2 + x3, random = ~ 1 | random, example.data)

partial.resid(y ~ x3, lme.model, example.data, return.data.frame = FALSE)

# Run mixed effects model using lmer()
library(lme4) 

lmer.model = lmer(y ~ x1 + x2 + x3 + (1 | random), example.data)

partial.resid(y ~ x3, lmer.model, example.data, return.data.frame = FALSE)

# Remove some values of x3 from data.frame
example.data[c(12, 23, 45), "x3"] = NA

# Run regular linear model using lm()
lm.model = lm(y ~ x1 + x2 + x3, example.data)

sum(!is.na(partial.resid(y ~ x3, lm.model, example.data)$x3)) # Should be 97

### SEM EXAMPLE ###

# Load model package
library(nlme)

# Load data from Shipley (2013)
data(shipley2013) 

shipley2013.modlist = list(
  
  lme(x2~x1, random = ~x1 | species, data = shipley2013),
  
  lme(x3~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x4~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x5~x3+x4, random = ~x3+x4 | species, data = shipley2013)
  
)

# Get partial residuals of x3 on x5 conditional on x4
resids.df = partial.resid(x5 ~ x3, shipley2013.modlist, shipley2013, list(lmeControl(opt = "optim")))

# Also returns raw residuals values for plotting in other packages
head(resids.df)
}
\name{partialResid}
\alias{partialResid}
\title{
Computing partial effects
}
\description{
Extracts partial residuals from a model or \code{psem} object for a given \code{x} and \code{y}.
}
\usage{
partialResid(formula., modelList, data = NULL)
}
\arguments{
  \item{formula.}{
  A formula where the \code{lhs} is the response and the \code{rhs} is the predictor whose partial effect is desired.
}
  \item{modelList}{
  A list of structural equations.
}
  \item{data}{
  A \code{data.frame} used to fit the equations.
}
}
\details{
This function computes the partial residuals of \code{y ~ x + Z} in a two-step procedure to remove the variation explained by \code{Z}: (1) remove \code{x} from the equation and model \code{y ~ Z}, and (2) replace \code{y} with \code{x} and model \code{x ~ Z}.
}
\value{
Returns a \code{data.frame} of residuals of \code{y ~ Z} called \code{yresids}, of \code{x ~ Z} called \code{xresids}.
}
\author{
Jon Lefcheck <jlefcheck@bigelow.org>
}
\seealso{
\code{\link{cerror}}
}
\examples{
# Generate data
dat <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))

# Build model
model <- lm(y ~ x1 + x2, dat)

# Compute partial residuals of y ~ x1
yresid <- resid(lm(y ~ x2, dat))

xresid <- resid(lm(x1 ~ x2, dat))

plot(yresid, xresid)

# Use partialResid
presid <- partialResid(y ~ x1, model)

plot(presid) # identical plot!
}

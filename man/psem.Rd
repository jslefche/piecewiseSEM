\name{psem}
\alias{psem}
\title{
Fitting piecewise structural equation models
}
\description{
Creates a list of structural equations.
}
\usage{
psem(..., data)
}
\arguments{
  \item{\dots}{
  A list of structural equations.
  }
  \item{data}{
  A \code{data.frame} used to fit the equations.
  }
}
\details{
\code{psem} takes a list of structural equations, which can be model objects of classes: \code{lm, glm, gls, pgls, sarlm, lme, glmmPQL, lmerMod, merModLmerTest, glmerMod}.

It also takes objects of class \code{formula, formula.cerror}, corresponding to additional variables to be included in the tests of directed separation (\code{X ~ 1}) or correlated errors (\code{X1 ~~ X2}).

The function optionally accepts data objects of classes: \code{matrix, data.frame, SpatialPointsDataFrame, comparative.data}, or these are derived internally from the structural equations.
}
\value{
Returns an object of class \code{psem}.
}
\author{
Jon Lefcheck <jlefcheck@bigelow.org>
}
\seealso{
\code{\link{\%~~\%}}
}
\examples{

}

\name{rsquared}
\alias{rsquared}
\title{
R-squared for linear regression
}
\description{
Returns (pseudo)-R2 values for all linear, general(ized) linear, and general(ized) linear mixed effects models.
}
\usage{
rsquared(modelList, method = NULL)
}
\arguments{
  \item{modelList}{
  A regression, or a list of structural equations.
}
  \item{method}{
  An optional vector of methods for use in calculating R2.
}
}
\details{


For mixed models, marginal R2 considers only the variance by the fixed effects, and the conditional R2 by both the fixed and random effects.
}
\value{
Returns a \code{data.frame} with the response, its family and link, and the R2 value. Mixed models also return marginal and conditional R2 values.
}
\references{
Nakagawa, Sinichi, Johnson, Paul C.D., and Holger Schielzeth. "The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisted and expanded." bioRxiv 095851 (2017).
}
\author{
Jon Lefcheck <jlefcheck@bigelow.org>
}
\examples{


}

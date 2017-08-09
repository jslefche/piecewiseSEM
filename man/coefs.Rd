\name{coefs}
\alias{coefs}
\title{
Extract path coefficients
}
\description{
Extracts (standardized) path coefficients from a \code{psem} object.
}
\usage{
coefs(modelList, intercepts = FALSE, standardize = TRUE)
}
\arguments{
  \item{modelList}{
  A list of structural equations.
}
  \item{intercepts}{
  Whether intercepts should be included in the coefficients table. Default is FALSE.
}
  \item{standardize}{
  Whether standardized path coefficients should be reported. Default is TRUE.    ~~Describe \code{standardize} here~~
}
}
\details{
P-values for models constructed using \code{\link{lme4::lme4}} are obtained using the Kenward-Roger approximation of the denominator degrees of freedom as implemented in the \code{\link{pbkrtest::pbkrtest}} package.
}
\value{
Returns a \code{data.frame} of coefficients, their standard errors, degrees of freedom, and significance tests.
}
\author{
Jon Lefcheck <jlefcheck@bigelow.org>
}
\seealso{
\code{\link{KRmodcomp}}
}
\examples{

}

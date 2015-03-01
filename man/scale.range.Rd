\name{scale.range}
\alias{scale.range}
\title{
Perform range standardization
}
\description{
Standardizes variables to interval [0-1].
}
\usage{
scale.range = function(x)
}
\arguments{
   \item{x}{
a numeric \code{vector}. 
}
}
\details{
A helper function for \code{get.sem.coefs} when \code{standardized = "range"}.
}
\value{
}
\references{
}
\author{
Jon Lefcheck
}
\note{
}
\seealso{
\code{\link{get.sem.coefs}}, \code{\link{scale}}
}
\examples{
x = runif(100, 0, 1)
scale.range(x)
}
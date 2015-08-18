\name{get.dag}

\alias{get.dag}

\title{
Generates adjacency matrix from list of structured equations
}

\description{
A helper function for \code{sem.missing.paths} that generates an adjacency matrix from a list of model formulae (generally recovered from a list of structured equations). 
}

\usage{
get.dag(formulaList)
}

\arguments{
  \item{formulaList}{
a list of model formula(e).
}
}

\details{
}

\value{
Returns a \code{matrix} of 0s and 1s, 0 indicating no relationship and 1 indicating a relationship, with columns being predictors and rows being responses.
}

\references{
}

\author{
Jon Lefcheck
}

\note{
}

\seealso{
}

\examples{
}
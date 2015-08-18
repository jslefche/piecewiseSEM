\name{get.basis.set}

\alias{get.basis.set}

\title{
Generates basis set for tests of direct separation
}

\description{
A helper function for \code{sem.missing.paths} that returns a list of (conditional) independence claims from an adjacency matrix generated from \code{get.dag}.
}

\usage{
get.basis.set(amat)
}

\arguments{
  \item{amat}{
an adjacency matrix.
}
}

\details{
}

\value{
Returns a list of independence claims. For each entry, the first argument is the variable whose partial regression slope should be zero, the second argument is the response, and the final entries are the variable(s) upon which the first argument is conditional.
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
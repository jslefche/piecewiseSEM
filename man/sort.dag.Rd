\name{sort.dag}

\alias{sort.dag}

\title{
Sorts adjacency matrix from parents to children
}

\description{
A helper function for \code{sem.basis.set} that takes the adjacency matrix from \code{get.dag} and orders it based on the appearance of variables in the SEM. It first puts variables that have no links (i.e., exogenous variables), then those that have only links to exogenous variables, then links to both exogenous and those that have links (i.e., endogenous), and finally to only endogenous.
}

\usage{
sort.dag(amat)
}

\arguments{
  \item{amat}{
adjacency matrix from \code{get.dag}
}
}

\value{
Returns a \code{matrix} of 0s and 1s, 0 indicating no relationship and 1 indicating a relationship, with columns being predictors and rows being responses.
}

\author{
Jon Lefcheck
}
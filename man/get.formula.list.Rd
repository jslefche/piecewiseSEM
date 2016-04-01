\name{get.formula.list}

\alias{get.formula.list}

\title{
Generates list of formula from a list of models.
}

\description{
A helper function for \code{sem.missing.paths} that retrieves a formula list from a list of structured equations.
} 

\usage{
get.formula.list(modelList)
}

\arguments{
  \item{modelList}{
a \code{list} of regressions representing the structural equation model.
}
}

\value{
Returns a \code{list} of objects of class \code{formula}.
}

\author{
Jon Lefcheck
}
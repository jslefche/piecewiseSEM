\name{get.random.formula}

\alias{get.random.formula}

\title{
Recovers random structure from a mixed model
}

\description{
A helper function for \code{sem.missing.paths} and \code{sem.partial.resid} that retrives (and modifies) the proper random structure for a mixed model to be passed to \code{update}.
}

\usage{
get.random.formula(model, rhs,  modelList, dropterms)
}

\arguments{
  \item{model}{
a linear model.
}
  \item{rhs}{
the fixed effects formula.
}
  \item{modelList}{
a \code{list} of regressions representing the structural equation model.
}
  \item{dropterms}{
terms to be dropped from the random structure (optional).
}
}

\value{
Returns a character vector corresponding to the updated random structure.
}
\author{
Jon Lefcheck
}

\seealso{
\code{\link[lme4]{findbars}}
}
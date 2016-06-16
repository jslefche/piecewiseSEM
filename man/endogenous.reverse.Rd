\name{endogenous.reverse}

\alias{endogenous.reverse}

\title{
Updates basis set under certain conditions.
}

\description{
A helper function for \code{sem.missing.paths} that reverses independence claims between intermediate endogenous variables fit to non-normal distributions.
} 

\usage{
endogenous.reverse(basis.set, modelList)
}

\arguments{
  \item{basis.set}{
a \code{list} of tests of conditional independence claims obtained from \code{sem.basis.set}.
}

\item{modelList}{
a \code{list} of regressions representing the structural equation model. 
}
}

\value{
Returns a \code{list} of conditional independence claims. The list contains duplicated entries, reversing the independence claims between endogenous variables that also predict other variables in the model *only if* those variables are fitted to a non-normal distribution.

This is a temporary workaround for GLM(M)s that do not adhere to the property of symmetricity of independence claims as a function of transformations via the link function. 
}

\author{
Jon Lefcheck
}

\examples{
\dontrun{
  data = data.frame(
    a = runif(50),
    b = runif(50, 0, 1),
    c = runif(50, 0, 1),
    d = runif(50)
  )
  
  modelList = list(
    glm(b ~ a, family = "binomial", data),
    glm(c ~ a, family = "binomial", data),
    lm(d ~ b + c, data)
  )
  
  (basis.set = sem.basis.set(modelList))
  
  endogenous.reverse(basis.set, modelList)
}
}
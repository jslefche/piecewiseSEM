\name{filter.exogenous}

\alias{filter.exogenous}

\title{
Filter exogenous variables from the basis set for SEM
}

\description{
Identifies exogenous variables (variables that have no paths leading to them) and removes them from the basis set when they appear as responses.
}

\usage{
filter.exogenous(modelList, basis.set, corr.errors, add.vars)
}

\arguments{

 \item{modelList}{
a \code{list} of regressions representing the structural equation model. 
}

  \item{basis.set}{
an \code{list} of vectors representing conditional independence statements.
}

  \item{corr.errors}{
a vector of variables with correlated errors (separated by "~~").
}

  \item{add.vars}{
a vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}

}

\details{
Helper function for \code{sem.basis.set} that removes independence claims between two exogenous (i.e., independent) variables, as the direction of causality cannot be reliably inferred from the SEM model list.
}

\value{
Returns a \code{list} of independence claims. Each entry in the \code{list} is a vector where the first entry is the predictor whose independence from the response is being evaluated, the second is the response, and remaining entries represent the variables on which the independence claim is conditional. E.g., 

  [1] \code{"x1" "y" "x2"} 

would be the missing path between \code{y ~ x1} conditional on \code{x2}. The basis set is constructed using the function \code{basiSet} in the \code{ggm} package.
}

\references{
Shipley, Bill. "A new inferential test for path models based on directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.
}

\author{
Jon Lefcheck
}

\examples{
# Load example data
data(shipley2009)

# Load model packages
library(lme4)
library(nlme)

# Create list of models 
shipley2009.modlist = list(
  
  lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Date ~ DD, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  lme(Growth ~ Date, random = ~1|site/tree, na.action = na.omit, 
      data = shipley2009),
  
  glmer(Live ~ Growth+(1|site)+(1|tree), 
        family=binomial(link = "logit"), data = shipley2009) 
  
)

# Filter exogenous variables
filter.exogenous(shipley2009.modlist)

# Compare with output from sem.basis.set
identical(
  
  sem.basis.set(shipley2009.modlist),
  
  filter.exogenous(shipley2009.modlist) 
  
  )
  
# Results are identical!

}
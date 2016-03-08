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
an optional \code{list} of vectors representing conditional independence statements.
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
A \code{list} of vectors representing conditional independence statements. The first string is the response, the second is the predictor, and the remaining strings are the variables (if any) upon which the indepedence claim are conditional.
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
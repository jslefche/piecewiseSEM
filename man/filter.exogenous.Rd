\name{filter.exogenous}
\alias{filter.exogenous}
\title{
Filter exogenous variables from the basis set
}
\description{
Identifies exogenous variables (variables that have no paths leading to them) and removes them from the basis set when they appear as responses.
}
\usage{
filter.exogenous(modelList, basis.set = NULL, corr.errors = NULL, 
add.vars = NULL)
}
\arguments{
 \item{modelList}{
A \code{list} of regressions representing the structural equation model. 
}
  \item{basis.set}{
A optional \code{list} of vectors representing conditional independence statements.
}
  \item{corr.errors}{
A vector of variables with correlated errors (separated by "~~").
}
  \item{add.vars}{
A vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
}
\details{

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
\note{
}
\examples{
data(Shipley) # Data from Shipley 2009

# Load model packages
library(lmerTest)
library(nlme)
# Create list of models 
Shipley.modlist = list(
  lme(DD~lat, random = ~1|site/tree, na.action = na.omit, 
  data = Shipley),
  lme(Date~DD, random = ~1|site/tree, na.action = na.omit, 
  data = Shipley),
  lme(Growth~Date, random = ~1|site/tree, na.action = na.omit, 
  data = Shipley),
  glmer(Live~Growth+(1|site)+(1|tree), 
  family=binomial(link = "logit"), data = Shipley) )

# Filter exogenous variables
filter.exogenous(Shipley.modlist)
# Compare with output from get.basis.set
identical(
  get.basis.set(Shipley.modlist),
  filter.exogenous(Shipley.modlist) )
# Results are identical (get.basis.set calls filter.exogenous)
}
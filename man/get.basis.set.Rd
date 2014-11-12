\name{get.basis.set}
\alias{get.basis.set}
\title{
Derive basis set of independence claims
}
\description{
Generates a list representing the (conditional) independence claims from a model list.
}
\usage{
get.basis.set(modelList, corr.errors = NULL, add.vars = NULL)
}
\arguments{
  \item{modelList}{
A \code{list} of regressions representing the structural equation model. 
}
  \item{corr.errors}{
A vector of variables with correlated errors (separated by "~~").
}
  \item{add.vars}{
A vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
}
\details{
Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{get.sem.coefs}.

The argument \code{add.vars} requires a vector of character strings corresponding to column names in the dataset used to construct the models in \code{modelList}. This is useful if comparing nested SEMs where one wishes to account for additional variables whose independence claims should be evaluated, but which do not have any hypothesized paths in the current SEM. The default is \code{NULL} (there is no additional independence claims that do not appear in the model list).
}
\value{
Returns a \code{list} of independence claims. Each entry in the \code{list} is a vector where the first entry is the predictor whose independence from the response is being evaluated, the second is the response, and remaining entries represent the variables on which the independence claim is conditional. 
}
\references{
Shipley, Bill. "A new inferential test for path models based on directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.
}
\author{
Jon Lefcheck
}
\note{
Unlike the functions \code{DAG} and \code{basiSet} in the \code{ggm} package, this function incorporates interactions if they included in the models.
}
\seealso{
\code{\link{DAG}}, \code{\link{basiSet}}
}
\examples{
# Fit model from Shipley (2009)
data(shipley2009) 

# Load model packages
library(lmerTest)
library(nlme)

# Create list of models 
shipley2009.modlist = list(
  lme(DD~lat, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  lme(Date~DD, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  lme(Growth~Date, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  glmer(Live~Growth+(1|site)+(1|tree), 
  family=binomial(link = "logit"), data = shipley2009) )

get.basis.set(shipley2009.modlist)

}
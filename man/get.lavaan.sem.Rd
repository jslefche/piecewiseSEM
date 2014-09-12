\name{get.lavaan.sem}
\alias{get.lavaan.sem}
\title{
Constructs variance-covariance based structural equation model.
}
\description{
Estimates variance-covariance based (traditional) SEM using \code{lavaan}.
}
\usage{
get.lavaan.sem(modelList, data, corr.errors = NULL, add.vars = NULL)
}
\arguments{
  \item{modelList}{
A list of regressions representing the structural equation model. 
}
  \item{data}{
A \code{data.frame} containing the full dataset for the SEM.
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
Returns a model object of class \code{lavaan}.
}
\references{
Grace, James B. Structural equation modeling and natural systems. Cambridge University Press, 2006.
}
\author{
Jon Lefcheck
}
\note{
}

\seealso{
\code{\link[lavaan]{sem}}
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
  
# Get variance-covariance based SEM
summary(get.lavaan.sem(Shipley.modlist, Shipley))
}
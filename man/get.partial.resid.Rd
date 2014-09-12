\name{get.partial.resid}
\alias{get.partial.resid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates partial residuals for two variables.
}
\description{
Extracts partial residuals for y ~ x | Z, where Z represents all other variables present in a structured equation upon which y and x are conditional.
}
\usage{
get.partial.resid(y, x, modelList)
}
\arguments{
   \item{y}{
The response variable whose partial relationship with x is desired.
}
  \item{x}{
The predictor variable whose partial relationship with y is desired. 
} 
 \item{modelList}{
A list of regressions representing the structural equation model. 
}
}
\details{

}
\value{
Returns a \code{data.frame} where the first column are the partial residuals of var1 ~ Z, and the second column is the partial residuals of var2 ~ Z.
}
\references{
Shipley, Bill. Cause and correlation in biology: a user's guide to path analysis, structural equations and causal inference. Cambridge University Press, 2002.
}
\author{
Jon Lefcheck
}
\note{
}
\seealso{

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
  
# Get path coefficients
get.partial.resid(y = "Date", x = "DD", Shipley.modlist) # Returns error
}
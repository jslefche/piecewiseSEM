\name{sem.lavaan}

\alias{sem.lavaan}

\title{
Constructs variance-covariance based SEM
}

\description{
Estimates variance-covariance based (traditional) structural equation model (SEM) using the \code{lavaan} package.
}

\usage{
sem.lavaan(modelList, data, corr.errors, add.vars, ...)
}

\arguments{

\item{modelList}{
a list of regressions representing the structural equation model. 
}

\item{data}{
a \code{data.frame} containing the full dataset for the SEM.
}

\item{corr.errors}{
a vector of variables with correlated errors (separated by "~~").
}

\item{add.vars}{
a vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}

\item{...}{
additional arguments passed to \code{lavaan}.
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
# Load model packages
library(lmerTest)
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
  
# Get variance-covariance based SEM
(lavaan.model = sem.lavaan(shipley2009.modlist, shipley2009))

summary(lavaan.model)

### 

# Load model packages
library(lmerTest)
library(nlme)

# Load data from Shipley (2013)
data(shipley2013) 

shipley2013.modlist = list(

  lme(x2~x1, random = ~x1 | species, data = shipley2013),
  
  lme(x3~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x4~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x5~x3+x4, random = ~x3+x4 | species, data = shipley2013)
  
  )
  
# Get variance-covariance based SEM
summary(sem.lavaan(shipley2013.modlist, shipley2013))
}
\name{piecewiseSEM-package}

\alias{piecewiseSEM-package}

\alias{piecewiseSEM}

\docType{package}

\title{
Piecewise Structural Equation Modeling
}

\description{
Implements piecewise structural equation modeling in R, complete with goodness-of-fit tests and retrieval of model coefficients. Compared with traditional variance-covariance based SEM, piecewise SEM allows for fitting of models to different distributions and/or incorporates hierarchical/nested random structures.

Supported model classes include: \code{lm}, \code{glm}, \code{glm.nb}, \code{gls}, \code{pgls}, \code{lme}, \code{glmmPQL}, and \code{merModLmerTest}.
}

\details{

\tabular{ll}{
Package: \tab piecewiseSEM\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2016-02-03\cr
Depends: \tab R (3.2.X), ggm, lavaan\cr
Suggests: \tab MASS, lmerTest, nlme\cr
License: \tab MIT\cr
}

The primary functions in the package are \code{\link{sem.fit}} which performs goodness-of-fit tests, and \code{\link{sem.coefs}} which returns path coefficients (standardized, if specified) and standard errors.
}

\author{
Jon Lefcheck <jslefche@vims.edu>
}

\references{
Shipley, Bill. "A new inferential test for path models based on directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.

Shipley, Bill. Cause and correlation in biology: a user's guide to path analysis, structural equations and causal inference. Cambridge University Press, 2002.

Shipley, Bill. "Confirmatory path analysis in a generalized multilevel context." Ecology 90.2 (2009): 363-368.

Shipley, Bill. "The AIC model selection method applied to path analytic models compared using a d-separation test." Ecology 94.3 (2013): 560-564.
}

\keyword{ package }
\seealso{
\code{\link{sem.fit}}
\code{\link{sem.coefs}}
}

\examples{
# Load example data
data(shipley2009)

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
  
# Get goodness-of-fit and AIC
sem.fit(shipley2009.modlist, shipley2009)

# Extract path coefficients
sem.coefs(shipley2009.modlist, shipley2009)
}
\name{sem.plot}

\alias{sem.plot}

\title{
Plotting of piecewise SEMs
}

\description{
Generates mandala plot for piecewise structural equation models.
}

\usage{
sem.fit(modelList, data, table, show.nonsig = TRUE, scaling = 10, ...)
}

\arguments{

\item{modelList}{
a \code{list} of regressions representing the structural equation model. 
}

\item{data}{
a \code{data.frame} used to construct the structured equations.
}

  \item{table}{
a table of model coefficients obtained from \code{sem.coefs}.
}

\item{show.nonsig}{
whether to show non-significant paths (P >= alpha). Default is \code{TRUE}.
}

\item{scaling}{
scaling coefficient for path thickness. Default is 10. If \code{NA} then all paths are equal width.
}

\item{alpha}{
significance threshold. Default is 0.05.
}

\item{...}{
additional arguments to \code{sem.coefs}.
}
}

\details{
Accepts either a list of structured equations and the data, or a table of coefficients obtained from \code{sem.coefs}. 
}

\value{
Returns a plot where variables are equidistantly spaced in a circle. Arrows indicate the directionality of the relationship. The arrow width indicates the strength of the effect (based on estimates from \code{sem.coefs}). Dashed and grey arrows indicate non-significant effects. Red arrows indicate negative effects.

}

\author{
Jon Lefcheck
}

\seealso{
\code{sem.coefs}
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

# Plot output
sem.plot(shipley2009.modlist, shipley2009)

### NOT RUN ###

# # Alternately get coefficient table and pass to sem.plot
# coef.table = sem.coefs(shipley2009.modlist, shipley2009)
#
#sem.plot(coef.table)

### END NOT RUN ###
}
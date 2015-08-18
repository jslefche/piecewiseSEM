\name{sem.fit}

\alias{sem.fit}

\title{
Goodness-of-fit tests for piecewise SEM
}

\description{
Tests independence claims and calculates Fisher's C statistic and associated p-value, and AIC and AICc, for a piecewise structural equation model (SEM).
}

\usage{
sem.fit(modelList, data, corr.errors, add.vars, grouping.vars, 
  top.level.vars, adjust.p, basis.set, pvalues.df, model.control, .progressBar) 
}

\arguments{

\item{modelList}{
a \code{list} of regressions representing the structural equation model. 
}

\item{data}{
a \code{data.frame} used to construct the structured equations.
}

  \item{corr.errors}{
a vector of variables with correlated errors (separated by "~~").
}

\item{add.vars}{
a vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}

\item{grouping.vars}{
an optional variable that represents the highest level of a hierarchical dataset.
}

\item{top.level.vars}{
an optional \code{vector} of variables that are identical at the highest level of the hierarchy.
}

\item{adjust.p}{
whether p-values degrees of freedom should be adjusted (see below). Default is \code{FALSE}.
}

\item{basis.set}{
provide an optional basis set.
}

\item{pvalues.df}{
an optional \code{data.frame} corresponding to p-values for independence claims.
}

\item{model.control}{
a \code{list} of model control arguments to be passed to d-sep models.
}

\item{.progressBar}{
enable optional text progress bar. Default is \code{TRUE}.
}

}

\details{
Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{sem.coefs}.

The argument \code{add.vars} requires a vector of character strings corresponding to column names in the dataset used to construct the models in \code{modelList}. This is useful if comparing nested SEMs where one wishes to account for additional variables whose independence claims should be evaluated, but which do not have any hypothesized paths in the current SEM. The default assumes there is no additional independence claims that do not appear in the model list.

For linear mixed effects models, p-values can be adjusted to accommodate the full model degrees of freedom using the argument \code{p.adjust = TRUE}. For more information, see Shipley 2013.

\code{top.level.vars} allows the user to optionally specify one or more top level variables for a hierarchical dataset, which will be summarized based on one or more grouping variables specified by \code{grouping.vars}. For example, consider a two-level hierarchy, where variables at the top level have identical values for each level of the grouping variable. If the response is a top level variable (is identically replicated for the grouping variable), this function takes a mean of the lower level variables for each level of the grouping variable, then runs the test of d-separation. If the response is fully replicated (occurs at the lower level), then no aggregation occurs.

The model controls in \code{model.control} will be assigned to any subsequent d-sep models derived from that model class. 

}

\value{
Returns a \code{list} with the following:

\item{missing.paths}{A \code{data.frame} where the first column is the independence claim, and the second through sixth columns the model estimates corresponding to the response variable in the independence claim. }

\item{fishers.c}{ A \code{data.frame} with the first entry corresponding to Fisher's C statistic, the second corresponding to the Chi-squared test degrees of freedom, and the third corresponding to the outcome (p-value) of the significance test derived from a Chi-squared distribution. }

\item{AIC}{ A \code{data.frame} where the first entry is the AIC score, and the second is the AICc score, and the third is the likelihood degrees of freedom (K). }

}

\references{
Shipley, Bill. "Confirmatory path analysis in a generalized multilevel context." Ecology 90.2 (2009): 363-368.

Shipley, Bill. "The AIC model selection method applied to path analytic models compared using a d-separation test." Ecology 94.3 (2013): 560-564.
}

\author{
Jon Lefcheck
}

\note{
}

\seealso{
\code{\link{ggm::dag}}, \code{\link{sem.missing.paths}}, \code{\link{sem.fisher.c}}, \code{\link{aic}}
}

\examples{

# Load model packages
library(lmerTest)
library(nlme)

# Load data from Shipley (2009)
data(shipley2009) 

# Create list of models 
shipley2009.modlist = list(

  lme(DD~lat, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  
  lme(Date~DD, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  
  lme(Growth~Date, random = ~1|site/tree, na.action = na.omit, 
  data = shipley2009),
  
  glmer(Live~Growth+(1|site)+(1|tree), 
  family=binomial(link = "logit"), data = shipley2009) 
  
  )
  
# Conduct goodness-of-fit tests
sem.fit(shipley2009.modlist, shipley2009)

### NOT RUN ###

# Repeat but use lme4 to construct all models

# shipley2009.modlist.lme4 = list(
#   
#   lmer(DD~lat + (1|site/tree), na.action = na.omit, 
#        data = shipley2009),
#   
#   lmer(Date~DD + (1|site/tree), na.action = na.omit, 
#        data = shipley2009),
#   
#   lmer(Growth~Date + (1|site/tree), na.action = na.omit, 
#        data = shipley2009),
#   
#   glmer(Live~Growth+(1|site)+(1|tree), 
#         family=binomial(link = "logit"), data = shipley2009) 
#   
# )

# sem.fit(shipley2009.modlist.lme4, shipley2009)

### NOT RUN ###

# Add new variable with no hypothesized links
set.seed(1)
add.var = rnorm(nrow(shipley2009), 50, 20)

# Test for independence
sem.fit(shipley2009.modlist, shipley2009, add.vars = c("add.var"))

###

# Add new variable grouped at site level
set.seed(2)
shipley2009$site.var = rep(rnorm(20, 10, 5), each = 95)

# Add a final model regressing at the level of the grouping variable
shipley2009.modlist2 = append(
  
  shipley2009.modlist, list(
  
  lm(site.var~Growth, data = aggregate(shipley2009, by = list(site = shipley2009[,"site"]), 
    mean, na.rm = TRUE))) 
  
  )

# Test for independence
sem.fit(shipley2009.modlist2, shipley2009, grouping.vars = c("site"), 
  top.level.vars = c("site.var"))

###

# Fit model from Shipley (2013)
data(shipley2013)

shipley2013.modlist = list(

  lme(x2~x1, random = ~x1 | species, data = shipley2013),
  
  lme(x3~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x4~x2, random = ~x2 | species, data = shipley2013),
  
  lme(x5~x3+x4, random = ~x3+x4 | species, data = shipley2013)
  
  )

## NOT RUN ##
# sem.fit(shipley2013.modlist, shipley2013) 
# Convergence error!
## NOT RUN ##

# Specify model controls by switching to old optimizer
sem.fit(shipley2013.modlist, shipley2013, model.control = list(lmeControl(opt = "optim")))

###

## NOT RUN ##
# Compare to output from ggm package
# # Load library
# library(ggm)
# 
# # Load data
# data(marks)
# # Generate direct acyclic graph (DAG)
# dag = DAG(mechanics ~ vectors+algebra, 
#           vectors ~ algebra,
#           statistics ~ algebra+analysis, 
#           analysis ~ algebra)
# # Run test of directed separation
# shipley.test(dag, cov(marks), n=88)
# 
# # Now create list of structured equations using lm()
# modelList = list(
#   lm(mechanics ~ vectors+algebra, data = marks),
#   lm( vectors ~ algebra, data = marks),
#   lm(statistics ~ algebra+analysis, data = marks),
#   lm(analysis ~ algebra, data = marks)
# )
# 
# # Test model fit
# sem.fisher.c(modelList, marks)
## NOT RUN ##

}
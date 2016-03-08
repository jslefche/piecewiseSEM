\name{sem.fit}

\alias{sem.fit}

\title{
Goodness-of-fit tests for piecewise SEM
}

\description{
Tests independence claims and calculates Fisher's C statistic and associated p-value, and AIC and AICc, for a piecewise structural equation model (SEM).
}

\usage{
sem.fit(modelList, data, conditional = FALSE, corr.errors = NULL, 
add.vars = NULL, grouping.vars = NULL, adjust.p = FALSE, basis.set, 
model.control = NULL, .progressBar = TRUE)
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

\item{adjust.p}{
whether p-values degrees of freedom should be adjusted (see below). Default is \code{FALSE}.
}

\item{basis.set}{
provide an optional basis set.
}

\item{conditional}{
whether conditional variables should be shown in the independence claim (unless the formula is fewer than 30 characters). Default is FALSE.
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
Independence claims are constructed based on how the variables are treated as in the model list. For example, if the indepedence claim includes a binary variable that is fit to a binomial distribution using an identity link, the function will evaluate the any claims using the same parameters. 

Similarly, for linear mixed effects models construted in \code{lme4} or \code{nlme}, varying slopes and intercepts are treated as in the model list. For example, if a variable is modeled with both a random slope and intercept in any model in the model list, that variable will be modeled with a random slope and intercept when evaluating all independence claims in which it appears. If slopes and intercepts vary for multiple variables, they will appear as such, even if they are conditional.

For models constructed with \code{lme4}, p-values are derived from a t-distribution with degrees of freedom approximated by the Kenward-Rogers correction, as implemented in the \code{pbkrtest} package.

For linear mixed effects models, p-values can be adjusted to accommodate the full model degrees of freedom using the argument \code{p.adjust = TRUE}. For more information, see Shipley 2013.

Variables with correlated errors have no direct relationship but rather are hypothesized to be driven by the same underlying factor. This covariance should be reflected as correlated errors (double-headed arrow). Correlated errors are specified using the syntax from the \code{lavaan} package: \code{var1 ~~ var2}. Variables with correlated errors are ignored in the basis set under the assumption that their correlations will be quantified later using the function \code{sem.coefs}.

The argument \code{add.vars} requires a vector of character strings corresponding to column names in the dataset used to construct the models in \code{modelList}. This is useful if comparing nested SEMs where one wishes to account for additional variables whose independence claims should be evaluated, but which do not have any hypothesized paths in the current SEM. The default assumes there is no additional independence claims that do not appear in the model list.

If the data is hierarchical dataset and one or more responses are identical for each level of a grouping factors -- artificially inflating the df -- users can summarize the dataset for each grouping factor(s) specified in the argument \code{grouping.vars}. For example, consider a two-level hierarchy, where variables at the top level have identical values for each level of the grouping variable. If the response is a top level variable (is identically replicated for the grouping variable), this function takes a mean of the lower level variables for each level of the grouping variable, then runs the test of d-separation. If the response is fully replicated (occurs at the lower level), then no aggregation occurs.
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
The model controls in \code{model.control} will be assigned to any subsequent d-sep models derived from that model class. 
}

\seealso{
\code{\link[ggm]{DAG}}, \code{\link{sem.missing.paths}}, \code{\link{sem.fisher.c}}, \code{\link{AIC}}, \code{\link{get_ddf_Lb}}
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

# Conduct goodness-of-fit tests
sem.fit(shipley2009.modlist, shipley2009)

### NOT RUN ###

# # Repeat but use lme4 to construct all models
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
#
# sem.fit(shipley2009.modlist.lme4, shipley2009)
#
# # Add new variable with no hypothesized links
# set.seed(1)
# add.var = rnorm(nrow(shipley2009), 50, 20)
#
# # Test for independence
# sem.fit(shipley2009.modlist, shipley2009, add.vars = c("add.var"))

###

# # Add new variable grouped at site level
# set.seed(2)
# shipley2009$site.var = rep(rnorm(20, 10, 5), each = 95)
# 
# # Add a final model regressing at the level of the grouping variable
# shipley2009.modlist2 = append(
#   
#   shipley2009.modlist, list(
#     
#     lm(site.var~Growth, data = aggregate(shipley2009, by = list(site = shipley2009[,"site"]), 
#                                          mean, na.rm = TRUE))) 
#   
# )
#
# # Test for independence
# sem.fit(shipley2009.modlist2, shipley2009, grouping.vars = c("site"))

###

# # Fit model from Shipley (2013)
# data(shipley2013)
# 
# shipley2013.modlist = list(
# 
#   lme(x2~x1, random = ~x1 | species, data = shipley2013),
# 
#   lme(x3~x2, random = ~x2 | species, data = shipley2013),
# 
#   lme(x4~x2, random = ~x2 | species, data = shipley2013),
# 
#   lme(x5~x3+x4, random = ~x3+x4 | species, data = shipley2013)
# 
# )
# 
# sem.fit(shipley2013.modlist, shipley2013) # Convergence error!
# 
# Specify model controls by switching to old optimizer
# sem.fit(shipley2013.modlist, shipley2013, model.control = list(lmeControl(opt = "optim")))

###

# # Compare to output from ggm package
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

### END NOT RUN ###
}
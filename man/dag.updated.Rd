\name{dag.updated}
\alias{dag.updated}
\title{
Directed acyclic graphs (DAGs)
}
\description{
Generates a matrix representing directed acylic graph (DAG) from a model list.
}
\usage{
dag.updated(modelList, add.vars = NULL)
}
\arguments{
  \item{modelList}{
A \code{list} of regressions representing the structural equation model. 
}
  \item{add.vars}{
A vector of additional variables whose independence claims should be evaluated, but which do not appear in the model list.
}
}
\details{
The argument \code{add.vars} requires a vector of character strings corresponding to column names in the dataset used to construct the models in \code{modelList}. This is useful if comparing nested SEMs where one wishes to account for additional variables whose independence claims should be evaluated, but which do not have any hypothesized paths in the current SEM. The default is \code{NULL} (there is no additional independence claims that do not appear in the model list).
}
\value{
Returns a adjacency matrix of the DAG, with rows and columns corresponding to the nodes of the graph, and a \code{1} if the variables are connected, or a \code{0} if they are not.
}
\references{
Shipley, Bill. "A new inferential test for path models based on directed acyclic graphs." Structural Equation Modeling 7.2 (2000): 206-218.
}
\author{
Jon Lefcheck
}
\note{
Unlike the function \code{DAG} in the \code{ggm} package, this function incorporates interactions if they included in the models.
}
\seealso{
\code{\link{DAG}} 
}
\examples{
data(shipley) # Data from Shipley 2009

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

dag.updated(Shipley.modlist)

}
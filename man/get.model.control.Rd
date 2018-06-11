\name{get.model.control}

\alias{get.model.control}

\title{
Formats model control arguments
}

\description{
A helper function for \code{sem.missing.paths} and \code{sem.partial.resid} that formats a \code{list} of model control arguments to pass to \code{update}.
}

\usage{
get.model.control(model, model.control)
}

\arguments{
  \item{model}{
a linear model.
}
  \item{model.control}{
a \code{list} of control parameters corresponding to the model classes found in the SEM.
}
}

\value{
Returns a \code{list} of model controls.
}


\author{
Jon Lefcheck
}

\seealso{
\code{\link{glm.control}}, \code{\link[nlme]{glsControl}}, \code{\link[nlme]{lmeControl}}, \code{\link[lme4]{lmerControl}}, \code{\link[lme4]{glmerControl}}
}
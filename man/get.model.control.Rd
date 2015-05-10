\name{get.model.control}

\alias{get.model.control}

\title{
Formats model control arguments
}

\description{
A helper function for \code{get.missing.paths} and \code{sem.partial.resid} that formats a \code{list} of model control arguments to pass to \code{update}.
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

\details{
}

\value{
Returns a \code{list} of model controls.
}

\references{
}

\author{
Jon Lefcheck
}

\note{
}

\seealso{
\code{\link{glm.control}}, \code{\link{glsControl}}, \code{\link{lmeControl}}, \code{\link{lmerControl}}, \code{\link{glmerControl}}
}

\examples{
}
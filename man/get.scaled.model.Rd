\name{get.scaled.model}

\alias{get.scaled.model}

\title{
Scales data frame
}

\description{
A helper function for \code{sem.coefs} that updates a model object with a scaled data frame (by mean and standard deviation, or range).
}

\usage{
get.scaled.model(model, newdata, modelList)
}

\arguments{

  \item{model}{
a linear model. 
}

\item{newdata}{
a \code{data.frame} from \code{get.scaled.data}.
}

  \item{modelList}{
a \code{list} of regressions representing the structural equation model. 
}
}

\value{
Returns a model object fit to the new data.
}

\author{
Jon Lefcheck
}

\seealso{
\code{\link{scale}}
}
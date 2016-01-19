\name{scale.data}

\alias{scale.data}

\title{
Scales data frame
}

\description{
A helper function for \code{sem.coefs} that provides a scaled data frame (by mean and standard deviation, or range).
}

\usage{
scale.data(modelList, data, standardize)
}

\arguments{

  \item{modelList}{
a list of regressions representing the structural equation model. 
}

\item{data}{
a \code{data.frame} used to construct the structured equations.
}

  \item{standardize}{
specifies the type of standardization to perform: \code{"none"}, \code{"scale"}, \code{"range"}.
}
}

\details{
Automatically calculates and scales variables that are transformed in the model formulae.
}

\value{
Returns a \code{data.frame} with the scaled (transformed) variables.
}

\references{
}

\author{
Jon Lefcheck
}

\note{
}

\seealso{
\code{\link{scale}}
}

\examples{
}
\name{get.scaled.data}

\alias{get.scaled.data}

\title{
Scales data frame
}

\description{
A helper function for \code{sem.coefs} that provides a scaled data frame (by mean and standard deviation, or range).
}

\usage{
get.scaled.data(modelList, data, standardize)
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

\author{
Jon Lefcheck
}

\seealso{
\code{\link{scale}}
}
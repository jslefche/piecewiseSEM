\name{coefs}
\alias{coefs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coefs(modelList, intercepts = FALSE, standardize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelList}{
%%     ~~Describe \code{modelList} here~~
}
  \item{intercepts}{
%%     ~~Describe \code{intercepts} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (modelList, intercepts = FALSE, standardize = TRUE) 
{
    if (!class(modelList) \%in\% c("psem", "list")) 
        modelList <- list(modelList)
    if (class(modelList) == "psem") 
        data <- modelList$data
    else data <- getData.(modelList)
    if (class(data) \%in\% c("SpatialPointsDataFrame")) 
        data <- data@data
    if (class(data) \%in\% c("comparative.data")) 
        data <- data$data
    modelList <- removeData(modelList, formulas = 2)
    if (standardize == TRUE) 
        ret <- stdCoefs(modelList, data, intercepts)
    else ret <- unstdCoefs(modelList, data, intercepts)
    ret <- cbind(ret, isSig(ret$P.Value))
    ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, 
        is.numeric))], 4)
    names(ret)[length(ret)] = ""
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{psem}
\alias{psem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
psem(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (...) 
{
    x <- list(...)
    idx <- which(sapply(x, function(y) all(class(y) \%in\% c("matrix", 
        "data.frame", "SpatialPointsDataFrame", "comparative.data"))))
    if (length(idx) > 0) {
        x <- x[c((1:length(x))[!1:length(x) \%in\% idx], idx)]
        names(x)[length(x)] <- "data"
    }
    else {
        x$data <- getData.(x)
    }
    if (any(is.na(names(x)))) {
        idx. <- which(is.na(names(x)))
        names(x)[idx.] <- idx.
    }
    evaluateClasses(x)
    formulaList <- listFormula(removeData(x, formulas = 1))
    if (any(duplicated(sapply(formulaList, function(y) all.vars.merMod(y)[1])))) 
        stop("Duplicate responses detected in the model list. Collapse into single multiple regression!", 
            call. = FALSE)
    class(x) <- "psem"
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

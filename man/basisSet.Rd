\name{basisSet}
\alias{basisSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
basisSet(modelList, direction = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelList}{
%%     ~~Describe \code{modelList} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
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
function (modelList, direction = NULL) 
{
    formulaList <- listFormula(modelList, remove = FALSE)
    amat <- Dag(formulaList)
    b <- lapply(1:nrow(amat), function(i) {
        lapply(i:ncol(amat), function(j) {
            if (amat[i, j] != 0 | i == j) 
                NULL
            else {
                cvars <- c(rownames(amat)[amat[, rownames(amat)[i], 
                  drop = FALSE] == 1], rownames(amat)[amat[, 
                  rownames(amat)[j], drop = FALSE] == 1])
                cvars <- cvars[!duplicated(cvars)]
                c(rownames(amat)[i], rownames(amat)[j], cvars)
            }
        })
    })
    b <- unlist(b, recursive = FALSE)
    b <- b[!sapply(b, is.null)]
    if (length(b) > 0) {
        b <- filterExisting(b, formulaList)
        b <- filterExogenous(b, amat)
        b <- removeCerror(b, formulaList)
        b <- reverseNonLin(modelList, b, amat, formulaList)
        if (!is.null(direction)) 
            b <- specifyDir(b, direction)
        b <- replaceTrans(modelList, b, amat)
    }
    return(b)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

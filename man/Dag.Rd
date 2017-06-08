\name{Dag}
\alias{Dag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Dag(formulaList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formulaList}{
%%     ~~Describe \code{formulaList} here~~
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
function (formulaList) 
{
    fList <- lapply(formulaList, function(i) if (any(class(i) \%in\% 
        "formula.cerror")) 
        NULL
    else i)
    fList <- fList[!sapply(fList, is.null)]
    fList <- lapply(fList, all.vars.merMod)
    vars <- unlist(fList)
    vars <- unname(vars[!duplicated(vars)])
    amat <- do.call(cbind, lapply(vars, function(i) {
        f <- which(sapply(fList, function(j) j[1] == i))
        if (length(f) == 0) 
            rep(0, length(vars))
        else {
            form <- fList[[which(sapply(fList, function(k) k[1] == 
                i))]]
            vars \%in\% form + 0
        }
    }))
    dimnames(amat) <- list(vars, vars)
    diag(amat) <- 0
    amat <- sortDag(amat, fList)
    if (cyclic(amat, fList)) 
        stop("Model is recursive. Remove feedback loops and re-run model")
    return(amat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

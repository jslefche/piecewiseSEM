\name{fisherC}
\alias{fisherC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fisherC(dTable, add.claims = NULL, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dTable}{
%%     ~~Describe \code{dTable} here~~
}
  \item{add.claims}{
%%     ~~Describe \code{add.claims} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
}
  \item{conserve}{
%%     ~~Describe \code{conserve} here~~
}
  \item{conditional}{
%%     ~~Describe \code{conditional} here~~
}
  \item{.progressBar}{
%%     ~~Describe \code{.progressBar} here~~
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
function (dTable, add.claims = NULL, direction = NULL, conserve = FALSE, 
    conditional = FALSE, .progressBar = FALSE) 
{
    if (class(dTable) == "psem") 
        dTable <- dSep(dTable, direction, conserve, conditional, 
            .progressBar)
    if (length(dTable) == 0) {
        Cstat <- 0
        DF <- 0
        P <- 1
    }
    else {
        ps <- dTable$P.Value
        if (!is.null(add.claims)) {
            ps <- ps + add.claims
            message("Fisher's C has been adjusted to include additional claims not shown in the tests of directed separation.")
        }
        if (any(ps == 0)) 
            ps <- ps + 1e-20
        Cstat <- -2 * sum(log(ps))
        DF <- 2 * length(ps)
        P <- 1 - pchisq(Cstat, DF)
    }
    ret <- data.frame(Fisher.C = Cstat, df = DF, P.Value = P)
    ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, 
        is.numeric))], 3)
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

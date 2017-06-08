\name{infCrit}
\alias{infCrit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
infCrit(modelList, Cstat = NULL, add.claims = NULL, direction = NULL, conserve = FALSE, conditional = FALSE, .progressBar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelList}{
%%     ~~Describe \code{modelList} here~~
}
  \item{Cstat}{
%%     ~~Describe \code{Cstat} here~~
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
function (modelList, Cstat = NULL, add.claims = NULL, direction = NULL, 
    conserve = FALSE, conditional = FALSE, .progressBar = FALSE) 
{
    if (is.null(Cstat)) 
        Cstat <- fisherC(modelList, add.claims, direction, conserve, 
            conditional, .progressBar)
    modelList <- removeData(modelList, formulas = 1)
    K <- do.call(sum, lapply(modelList, function(i) attr(logLik(i), 
        "df")))
    n.obs <- min(sapply(modelList, nobs.))
    pwAIC <- as.numeric(Cstat[1] + 2 * K)
    pwAICc <- pwAIC * (n.obs/(n.obs - K - 1))
    pwBIC <- as.numeric(Cstat[1] + log(n.obs) * K)
    ret <- data.frame(AIC = pwAIC, AICc = pwAICc, BIC = pwBIC, 
        K = K, n = n.obs)
    ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, 
        is.numeric))], 3)
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
